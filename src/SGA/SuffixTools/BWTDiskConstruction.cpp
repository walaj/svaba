//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BWTDiskConstruction - Implementation of the
// BWT disk construction algorithm of Ferragina,
// Gagie, Manzini. See Lightweight Data Indexing
// and Compression in External Memory
//
#include "BWTDiskConstruction.h"
#include "Util.h"
#include "BWTWriter.h"
#include "BWTReader.h"
#include "SAWriter.h"
#include "SAReader.h"
#include "GapArray.h"
#include "RankProcess.h"
#include "SequenceProcessFramework.h"
#include "BWTCABauerCoxRosone.h"
#include <unistd.h>

// Definitions and structures
static const bool USE_GZ = false;
static const int BWT_SAMPLE_RATE = 512;

struct MergeItem
{
    int64_t start_index;
    int64_t end_index;
    std::string reads_filename;
    std::string bwt_filename;
    std::string sai_filename;

    friend std::ostream& operator<<(std::ostream& out, const MergeItem& item)
    {
        out << "[" << item.start_index << "," << item.end_index << "] " << item.bwt_filename;
        out << " " << item.sai_filename << " " << item.reads_filename;
        return out;
    }
};
typedef std::vector<MergeItem> MergeVector;

// Function declarations
int64_t merge(SeqReader* pReader, 
              const MergeItem& item1, const MergeItem& item2, 
              const std::string& bwt_outname, const std::string& sai_outname,
              bool doReverse, int numThreads, int storageLevel);

// Initial BWT construction algorithms
MergeVector computeInitialSAIS(const BWTDiskParameters& parameters); 
MergeVector computeInitialBCR(const BWTDiskParameters& parameters); 

//
void writeMergedIndex(const BWT* pBWTInternal, const MergeItem& externalItem, 
                      const MergeItem& internalItem, const std::string& bwt_outname,
                      const std::string& sai_outname, const GapArray* pGapArray);

void writeRemovalIndex(const BWT* pBWTInternal, const std::string& sai_inname,
                       const std::string& bwt_outname, const std::string& sai_outname, 
                       size_t num_strings_remove, size_t num_symbols_remove,
                       const GapArray* pGapArray);

void computeGapArray(SeqReader* pReader, size_t n, const BWT* pBWT, bool doReverse, 
                     int numThreads, GapArray* pGapArray, bool removeMode,
                     size_t& num_strings_read, size_t& num_symbols_read);

//
std::string makeTempName(const std::string& prefix, int id, const std::string& extension);
std::string makeFilename(const std::string& prefix, const std::string& extension);


// The algorithm is as follows. We create M BWTs for subsets of 
// the input reads. These are created independently and written
// to disk. They are then merged either sequentially or pairwise
// to create the final BWT
void buildBWTDisk(const BWTDiskParameters& parameters)
{
    // Build the initial bwts for subsets of the data
    MergeVector mergeVector;
    if(parameters.bUseBCR)
        mergeVector = computeInitialBCR(parameters);
    else
        mergeVector = computeInitialSAIS(parameters);

    // Phase 2: Pairwise merge the BWTs
    int groupID = mergeVector.size(); // Initial the name of the next intermediate bwt
    int round = 1;
    MergeVector nextMergeRound;
    while(mergeVector.size() > 1)
    {
        std::cout << "Starting round " << round << "\n";
        SeqReader* pReader = new SeqReader(parameters.inFile);
        SeqRecord record;

        for(size_t i = 0; i < mergeVector.size(); i+=2)
        {
            if(i + 1 != mergeVector.size())
            {
                std::string bwt_merged_name = makeTempName(parameters.outPrefix, groupID, parameters.bwtExtension);
                std::string sai_merged_name = makeTempName(parameters.outPrefix, groupID, parameters.saiExtension);

                MergeItem item1 = mergeVector[i];
                MergeItem item2 = mergeVector[i+1];

                // Perform the actual merge
                int64_t curr_idx = merge(pReader, item1, item2, 
                                         bwt_merged_name, sai_merged_name, 
                                         parameters.bBuildReverse, parameters.numThreads, parameters.storageLevel);

                // pReader now points to the end of item1's block of 
                // reads. Skip item2's reads
                assert(curr_idx == item2.start_index);
                while(curr_idx <= item2.end_index)
                {
                    bool eof = !pReader->get(record);
                    assert(!eof);
                    (void)eof;
                    ++curr_idx;
                }

                // Create the merged mergeItem to use in the next round
                MergeItem merged;
                merged.start_index = item1.start_index;
                merged.end_index = item2.end_index;
                merged.bwt_filename = bwt_merged_name;
                merged.sai_filename = sai_merged_name;
                nextMergeRound.push_back(merged);

                // Done with the temp files, remove them
                unlink(item1.bwt_filename.c_str());
                unlink(item2.bwt_filename.c_str());
                unlink(item1.sai_filename.c_str());
                unlink(item2.sai_filename.c_str());

                ++groupID;
            }
            else
            {
                // Singleton, pass through to the next round
                nextMergeRound.push_back(mergeVector[i]);
            }
        }
        delete pReader;
        mergeVector.clear();
        mergeVector.swap(nextMergeRound);
        ++round;
    }
    assert(mergeVector.size() == 1);

    // Done, rename the files to their final name
    std::stringstream bwt_ss;
    bwt_ss << parameters.outPrefix << parameters.bwtExtension << (USE_GZ ? ".gz" : "");
    std::string bwt_final_filename = bwt_ss.str();
    rename(mergeVector.front().bwt_filename.c_str(), bwt_final_filename.c_str());

    std::stringstream sai_ss;
    sai_ss << parameters.outPrefix << parameters.saiExtension << (USE_GZ ? ".gz" : "");
    std::string sai_final_filename = sai_ss.str();
    rename(mergeVector.front().sai_filename.c_str(), sai_final_filename.c_str());
}

// Compute the initial BWTs for the input file split into blocks of records using the SAIS algorithm
MergeVector computeInitialSAIS(const BWTDiskParameters& parameters)
{
    SeqReader* pReader = new SeqReader(parameters.inFile);
    SeqRecord record;

    int groupID = 0;
    size_t numReadTotal = 0;

    MergeVector mergeVector;
    MergeItem mergeItem;
    mergeItem.start_index = 0;

    // Phase 1: Compute the initial BWTs
    ReadTable* pCurrRT = new ReadTable;
    bool done = false;
    while(!done)
    {
        done = !pReader->get(record);

        if(!done)
        {
            // the read is valid
            SeqItem item = record.toSeqItem();
            if(parameters.bBuildReverse)
                item.seq.reverse();
            pCurrRT->addRead(item);
            ++numReadTotal;
        }

        if(pCurrRT->getCount() >= parameters.numReadsPerBatch || (done && pCurrRT->getCount() > 0))
        {
            // Compute the SA and BWT for this group
            SuffixArray* pSA = new SuffixArray(pCurrRT, 1);

            // Write the BWT to disk                
            std::string bwt_temp_filename = makeTempName(parameters.outPrefix, groupID, parameters.bwtExtension);
            pSA->writeBWT(bwt_temp_filename, pCurrRT);

            std::string sai_temp_filename = makeTempName(parameters.outPrefix, groupID, parameters.saiExtension);
            pSA->writeIndex(sai_temp_filename);

            // Push the merge info
            mergeItem.end_index = numReadTotal - 1; // inclusive
            mergeItem.reads_filename = parameters.inFile;
            mergeItem.bwt_filename = bwt_temp_filename;
            mergeItem.sai_filename = sai_temp_filename;
            mergeVector.push_back(mergeItem);

            // Cleanup
            delete pSA;

            // Start the new group
            mergeItem.start_index = numReadTotal;
            ++groupID;
            pCurrRT->clear();
        }
    }
    delete pCurrRT;
    delete pReader;
    return mergeVector;
}

// Compute the initial BWTs for the input file split into blocks of records using the BCR algorithm
MergeVector computeInitialBCR(const BWTDiskParameters& parameters)
{
    SeqReader* pReader = new SeqReader(parameters.inFile);
    SeqRecord record;

    int groupID = 0;
    size_t numReadTotal = 0;

    MergeVector mergeVector;
    MergeItem mergeItem;
    mergeItem.start_index = 0;

    // Phase 1: Compute the initial BWTs
    DNAEncodedStringVector readSequences;
    bool done = false;
    while(!done)
    {
        done = !pReader->get(record);

        if(!done)
        {
            // the read is valid
            SeqItem item = record.toSeqItem();
            if(parameters.bBuildReverse)
                item.seq.reverse();
            readSequences.push_back(item.seq.toString());
            ++numReadTotal;
        }

        if(readSequences.size() >= parameters.numReadsPerBatch || (done && readSequences.size() > 0))
        {
            std::string bwt_temp_filename = makeTempName(parameters.outPrefix, groupID, parameters.bwtExtension);
            std::string sai_temp_filename = makeTempName(parameters.outPrefix, groupID, parameters.saiExtension);
            BWTCA::runBauerCoxRosone(&readSequences, bwt_temp_filename, sai_temp_filename);

            // Push the merge info
            mergeItem.end_index = numReadTotal - 1; // inclusive
            mergeItem.reads_filename = parameters.inFile;
            mergeItem.bwt_filename = bwt_temp_filename;
            mergeItem.sai_filename = sai_temp_filename;
            mergeVector.push_back(mergeItem);

            // Start the new group
            mergeItem.start_index = numReadTotal;
            ++groupID;
            readSequences.clear();
        }
    }
    delete pReader;
    return mergeVector;
}

// Merge the indices for the two independent sets of reads in readsFile1 and readsFile2
void mergeIndependentIndices(const std::string& readsFile1, const std::string& readsFile2, 
                             const std::string& outPrefix, const std::string& bwt_extension, 
                             const std::string& sai_extension, bool doReverse, int numThreads, int storageLevel)
{
    MergeItem item1;
    std::string prefix1 = stripGzippedExtension(readsFile1);
    item1.reads_filename = readsFile1;
    item1.bwt_filename = makeFilename(prefix1, bwt_extension);
    item1.sai_filename = makeFilename(prefix1, sai_extension);
    item1.start_index = 0;
    item1.end_index = -1; // this tells merge to read the entire file

    MergeItem item2;
    std::string prefix2 = stripGzippedExtension(readsFile2);
    item2.reads_filename = readsFile2;
    item2.bwt_filename = makeFilename(prefix2, bwt_extension);
    item2.sai_filename = makeFilename(prefix2, sai_extension);
    item2.start_index = 0;
    item2.end_index = -1;

    // Prepare the input filehandle
    SeqReader* pReader = new SeqReader(item1.reads_filename);
    
    // Build the outnames
    std::string bwt_merged_name = makeFilename(outPrefix, bwt_extension);
    std::string sai_merged_name = makeFilename(outPrefix, sai_extension);

    // Perform the actual merge
    merge(pReader, item1, item2, bwt_merged_name, sai_merged_name, doReverse, numThreads, storageLevel);
    delete pReader;
}

// Construct new indices without the reads in readsToRemove
void removeReadsFromIndices(const std::string& allReadsPrefix, const std::string& readsToRemove,
                             const std::string& outPrefix, const std::string& bwt_extension, 
                             const std::string& sai_extension, bool doReverse, int numThreads)
{
    std::string bwt_filename = makeFilename(allReadsPrefix, bwt_extension);
    std::string sai_filename = makeFilename(allReadsPrefix, sai_extension);

    // Prepare the input filehandle
    SeqReader* pReader = new SeqReader(readsToRemove);
    
    // Build the outnames
    std::string bwt_out_name = makeFilename(outPrefix, bwt_extension);
    std::string sai_out_name = makeFilename(outPrefix, sai_extension);

    // Compute the gap array
    BWT* pBWT = new BWT(bwt_filename, BWT_SAMPLE_RATE);

    // Boolean gap array
    GapArray* pGapArray = createGapArray(1);

    size_t num_strings_remove;
    size_t num_symbols_remove;
    computeGapArray(pReader, (size_t)-1, pBWT, doReverse, numThreads, pGapArray, true, num_strings_remove, num_symbols_remove);

    //writeRemovalIndex();
    writeRemovalIndex(pBWT, sai_filename, bwt_out_name, sai_out_name, num_strings_remove, num_symbols_remove, pGapArray);

    // Perform the actual merge
    //merge(pReader, item1, item2, bwt_merged_name, sai_merged_name, doReverse, numThreads);

    delete pGapArray;
    delete pReader;
    delete pBWT;
}

// Merge two readsFiles together
void mergeReadFiles(const std::string& readsFile1, const std::string& readsFile2, const std::string& outPrefix)
{
    // If the outfile is the empty string, append the reads in readsFile2 into readsFile1
    // otherwise cat the files together
    std::ostream* pWriter;
    if(outPrefix.empty())
    {
        pWriter = createWriter(readsFile1, std::ios_base::out | std::ios_base::app);
    }
    else
    {
        bool both_fastq = isFastq(readsFile1) && isFastq(readsFile2);
        bool both_gzip = isGzip(readsFile1) && isGzip(readsFile2);
        std::string extension = both_fastq ? ".fastq" : ".fa";
        if(both_gzip)
            extension.append(".gz");
        std::string out_filename = outPrefix + extension;
        pWriter = createWriter(out_filename);

        // Copy reads1 to the outfile
        SeqReader reader(readsFile1);
        SeqRecord record;
        while(reader.get(record))
            record.write(*pWriter);
    }

    // Copy reads2 to writer
    SeqReader reader(readsFile2);
    SeqRecord record;
    while(reader.get(record))
        record.write(*pWriter);
    delete pWriter;
}

// Compute the gap array for the first n items in pReader
void computeGapArray(SeqReader* pReader, size_t n, const BWT* pBWT, bool doReverse, int numThreads, GapArray* pGapArray, 
                     bool removeMode, size_t& num_strings_read, size_t& num_symbols_read)
{
    // Create the gap array
    size_t gap_array_size = pBWT->getBWLen() + 1;
    pGapArray->resize(gap_array_size);

    // The rank processor calculates the rank of every suffix of a given sequence
    // and returns a vector of ranks. The postprocessor takes in the vector
    // and updates the gap array
    RankPostProcess postProcessor(pGapArray);
    size_t numProcessed = 0;
    if(numThreads <= 1)
    {
        RankProcess processor(pBWT, pGapArray, doReverse, removeMode);

        numProcessed = 
           SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
                                                            RankResult, 
                                                            RankProcess, 
                                                            RankPostProcess>(*pReader, &processor, &postProcessor, n);
    }
    else
    {
        typedef std::vector<RankProcess*> RankProcessVector;
        RankProcessVector rankProcVec;
        for(int i = 0; i < numThreads; ++i)
        {
            RankProcess* pProcess = new RankProcess(pBWT, pGapArray, doReverse, removeMode);
            rankProcVec.push_back(pProcess);
        }
    
        numProcessed = 
           SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
                                                              RankResult, 
                                                              RankProcess, 
                                                              RankPostProcess>(*pReader, rankProcVec, &postProcessor, n);

        for(int i = 0; i < numThreads; ++i)
            delete rankProcVec[i];
    }

    num_strings_read = postProcessor.getNumStringsProcessed();
    num_symbols_read = postProcessor.getNumSymbolsProcessed();
    assert(n == (size_t)-1 || (numProcessed == n));
}

// Merge a pair of BWTs using disk storage
// Precondition: pReader is positioned at the start of the read block for item1
int64_t merge(SeqReader* pReader,
              const MergeItem& item1, const MergeItem& item2,
              const std::string& bwt_outname, const std::string& sai_outname,
              bool doReverse, int numThreads, int storageLevel)
{
    std::cout << "Merge1: " << item1 << "\n";
    std::cout << "Merge2: " << item2 << "\n";

    // Load the bwt of item2 into memory as the internal bwt
    BWT* pBWTInternal = new BWT(item2.bwt_filename, BWT_SAMPLE_RATE);
    
    // If end_index is -1, calculate the ranks for every sequence in the file
    // otherwise only calculate the rank for the next (end_index - start_index + 1) sequences
    size_t n = (item1.end_index == -1) ? -1 : item1.end_index - item1.start_index + 1;
    
    // Calculate the rank of every read from item1.start_index to item1.end_index
    // and increment the gap counts
    int64_t curr_idx = item1.start_index;
    
    // Compute the gap/rank array
    GapArray* pGapArray = createGapArray(storageLevel);
    size_t num_strings_read = 0;
    size_t num_symbols_read = 0;
    computeGapArray(pReader, n, pBWTInternal, doReverse, numThreads, pGapArray, 
                    false, num_strings_read, num_symbols_read);

    assert(n == (size_t)-1 || (num_strings_read == n));

    // At this point, the gap array has been calculated for all the sequences
    curr_idx += num_strings_read;
    assert(item1.end_index == -1 || (curr_idx == item1.end_index + 1 && curr_idx == item2.start_index));

    // Write the merged BWT/SAI to disk
    writeMergedIndex(pBWTInternal, item1, item2, bwt_outname, sai_outname, pGapArray);

    delete pBWTInternal;
    delete pGapArray;
    return curr_idx;
}

// Merge the internal and external BWTs and the SAIs
void writeMergedIndex(const BWT* pBWTInternal, const MergeItem& externalItem, 
                      const MergeItem& internalItem, const std::string& bwt_outname,
                      const std::string& sai_outname, const GapArray* pGapArray)
{
    IBWTWriter* pBWTWriter = BWTWriter::createWriter(bwt_outname);
    IBWTReader* pBWTExtReader = BWTReader::createReader(externalItem.bwt_filename);
    
    SAWriter saiWriter(sai_outname);
    SAReader saiExtReader(externalItem.sai_filename);
    SAReader saiIntReader(internalItem.sai_filename);

    // Calculate and write header values
    size_t disk_strings;
    size_t disk_symbols;
    BWFlag flag;
    pBWTExtReader->readHeader(disk_strings, disk_symbols, flag);

    size_t total_strings = disk_strings + pBWTInternal->getNumStrings();
    size_t total_symbols = disk_symbols + pBWTInternal->getBWLen();
    pBWTWriter->writeHeader(total_strings, total_symbols, BWF_NOFMI);
    
    // Discard the first header each sai
    size_t discard1, discard2;
    saiExtReader.readHeader(discard1, discard2);
    saiIntReader.readHeader(discard1, discard2);

    // Write the header of the SAI which is just the number of strings and elements in the SAI
    saiWriter.writeHeader(total_strings, total_strings);

    // Calculate and write the actual string
    // The semantics of the gap array are that we need to write gap_array[i]
    // symbols to the stream before writing bwtInternal[i]
    // Each time a '$' symbol is read, it signals the end of some read. We
    // output one element of the sai from corresponding internal or external
    // sai file.
    size_t num_bwt_wrote = 0;
    size_t num_sai_wrote = 0;
    for(size_t i = 0; i < pGapArray->size(); ++i)
    {
        size_t v = pGapArray->get(i);
        for(size_t j = 0; j < v; ++j)
        {
            char b = pBWTExtReader->readBWChar();
            assert(b != '\n');
            pBWTWriter->writeBWChar(b);
            ++num_bwt_wrote;
            
            if(b == '$')
            {
                // The external indices only need to be copied
                SAElem e = saiExtReader.readElem(); 
                saiWriter.writeElem(e);
                ++num_sai_wrote;
            }
        }
        
        // If this is the last entry in the gap array, do not output a symbol from
        // the internal BWT
        if(i != pBWTInternal->getBWLen())
        {
            char b = pBWTInternal->getChar(i);
            pBWTWriter->writeBWChar(b);
            ++num_bwt_wrote;

            if(b == '$')
            {
                // The internal indices need to be offset
                // by the number of strings in the external collection
                SAElem e = saiIntReader.readElem(); 

                uint64_t id = e.getID();
                id += disk_strings;
                e.setID(id);
                
                saiWriter.writeElem(e);
                ++num_sai_wrote;
            }
        }
    }
    if(num_bwt_wrote != total_symbols)
    {
        printf("Error expected to write %zu symbols, actually wrote %zu\n", total_symbols, num_bwt_wrote);
        assert(num_bwt_wrote == total_symbols);
    }
        
    assert(num_sai_wrote == total_strings);
    
    // Ensure we read the entire bw string from disk
    char last = pBWTExtReader->readBWChar();
    assert(last == '\n');
    (void)last;

    // Finalize the BWT disk file
    pBWTWriter->finalize();

    delete pBWTExtReader;
    delete pBWTWriter;
}

// Write a new BWT and SAI that skips the elements marked
// by the gap array. This is used to remove entire strings from the 
// index
void writeRemovalIndex(const BWT* pBWTInternal, const std::string& sai_inname,
                       const std::string& bwt_outname, const std::string& sai_outname, 
                       size_t num_strings_remove, size_t num_symbols_remove,
                       const GapArray* pGapArray)
{
    IBWTWriter* pBWTWriter = BWTWriter::createWriter(bwt_outname);
    
    SAWriter* pSAIWriter = new SAWriter(sai_outname);
    SAReader* pSAIReader = new SAReader(sai_inname);

    // Calculate and write header values
    assert(num_strings_remove <= pBWTInternal->getNumStrings());
    assert(num_symbols_remove <= pBWTInternal->getBWLen());
    size_t input_strings = pBWTInternal->getNumStrings();
    size_t input_symbols = pBWTInternal->getBWLen();

    size_t output_strings = input_strings - num_strings_remove;
    size_t output_symbols = input_symbols - num_symbols_remove;

    pBWTWriter->writeHeader(output_strings, output_symbols, BWF_NOFMI);
    
    // Discard the header of the sai
    size_t discard1, discard2;
    pSAIReader->readHeader(discard1, discard2);

    // Write the header of the SAI which is just the number of strings and elements in the SAI
    pSAIWriter->writeHeader(output_strings, output_strings);

    // We need to fix the IDs in the SAI. We do this by tracking
    // the number of IDs that are removed that are lower than a given
    // id.
    std::vector<int> id_vector(pBWTInternal->getNumStrings(), 0);

    // Calculate and write the actual string
    // The gap array marks the symbols that should be removed. We
    // iterate over the gap array, if the value is zero, we output 
    // the BWT symbol. Otherwise we skip the number of elements
    // indicated by the gap. 
    size_t num_bwt_wrote = 0;
    size_t num_sai_wrote = 0;
    size_t i = 0;
    while(i < input_symbols)
    {
        size_t v = pGapArray->get(i);
        // If v is zero we output a single symbol, otherwise
        // we skip v elements. 
        size_t num_to_read = (v == 0) ? 1 : v;
        for(size_t j = 0; j < num_to_read; ++j)
        {
            char b = pBWTInternal->getChar(i);
            if(b == '$')
            {
                SAElem e = pSAIReader->readElem(); 
            
                // This element of the SAI will be removed,
                // increment the id vector so we can correct
                // the indices later
                if(v > 0)
                {
                    id_vector[e.getID()]++;
                }
            }
            
            // output
            if(v == 0)
            {
                pBWTWriter->writeBWChar(b);
                assert(b != '\n');
                ++num_bwt_wrote;
            }
        }
        i += num_to_read;
    }
    
    if(num_bwt_wrote != output_symbols)
    {
        printf("Error expected to write %zu symbols, actually wrote %zu\n", output_symbols, num_bwt_wrote);
        assert(num_bwt_wrote == output_symbols);
    }

    // Finalize the BWT disk file
    pBWTWriter->finalize();

    // Accumulate the counts for each id
    int prev = 0;
    for(size_t i = 0; i < id_vector.size(); ++i)
    {
        id_vector[i] += prev;
        prev = id_vector[i];
    }

    // Read the SAI file again, writing out the
    // corrected IDs that are not marked for removal
    delete pSAIReader;
    pSAIReader = new SAReader(sai_inname);
    pSAIReader->readHeader(discard1, discard2);
    for(size_t i = 0; i < input_strings; ++i)
    {
        SAElem e = pSAIReader->readElem();
        uint64_t id = e.getID();
        int prev_count = (id > 0) ? id_vector[id - 1] : 0;
        int count = id_vector[id];
        if(count == prev_count)
        {
            // the current element should be output, it wasn't marked for removal
            id -= count;
            e.setID(id);
            pSAIWriter->writeElem(e);
            ++num_sai_wrote;
        }
    }
    assert(num_sai_wrote == output_strings);

    delete pSAIReader;
    delete pSAIWriter;
    delete pBWTWriter;
}

//
std::string makeTempName(const std::string& prefix, int id, const std::string& extension)
{
    std::stringstream ss;
    ss << prefix << ".temp-" << id << extension << (USE_GZ ? ".gz" : "");
    return ss.str();
}

//
std::string makeFilename(const std::string& prefix, const std::string& extension)
{
    std::stringstream ss;
    ss << prefix << extension << (USE_GZ ? ".gz" : "");
    return ss.str();
}

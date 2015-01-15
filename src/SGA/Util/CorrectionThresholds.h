//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// CorrectionThresholds.h - Singeleton class containing
// cutoff values for determining whether a 
// base should be corrected or not
//
#ifndef CORRECTION_THRESHOLDS_H
#define CORRECTION_THRESHOLDS_H

class CorrectionThresholds
{
    public:
        CorrectionThresholds();
        static CorrectionThresholds& Instance();
        
        // Set the base minimum support level (for high-quality reads)
        void setBaseMinSupport(int ms);

        int getMinSupportHighQuality() { return m_minSupportHighQuality; }
        int getMinSupportLowQuality() { return m_minSupportLowQuality; }
        int getHighQualityCutoff() { return m_highQualityCutoff; }

        // Returns the support required for a base with phred score phred
        int getRequiredSupport(int phred);

    private:
        int m_highQualityCutoff;
        int m_minSupportLowQuality;
        int m_minSupportHighQuality;
};

#endif

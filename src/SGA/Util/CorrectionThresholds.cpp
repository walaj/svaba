//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// CorrectionThresholds - Cutoff values for determining
// whether a particular base should be corrected or not
//
#include "CorrectionThresholds.h" 

CorrectionThresholds::CorrectionThresholds()
{
    // Set defaults
    m_minSupportLowQuality = 4;
    m_minSupportHighQuality = 3;
    m_highQualityCutoff = 20;
}

CorrectionThresholds& CorrectionThresholds::Instance()
{
    static CorrectionThresholds instance;
    return instance;
}

void CorrectionThresholds::setBaseMinSupport(int ms)
{
    m_minSupportHighQuality = ms;
    m_minSupportLowQuality = ms + 1;
}

//
int CorrectionThresholds::getRequiredSupport(int phred)
{
    int threshold = m_minSupportLowQuality;
    if(phred >= m_highQualityCutoff)
        threshold = m_minSupportHighQuality;
    return threshold;
}


"""
    data for connecting src_tags with their corresponding reference information
"""


# direct mapping of src_tag to the full reference info
src_tag_to_reference = {
    'zhou1016': 'Zhou, Z., Shen, X., Tu, J. & Zhu, Z.-J. Large-Scale Prediction of Collision Cross-Section Values for Metabolites in Ion Mobility-Mass Spectrometry. Anal. Chem. 88, 11084–11091 (2016).',
    'zhou0817': 'Zhou, Z., Tu, J., Xiong, X., Shen, X. & Zhu, Z.-J. LipidCCS: Prediction of Collision Cross-Section Values for Lipids with High Precision To Support Ion Mobility–Mass Spectrometry-Based Lipidomics. Anal. Chem. 89, 9559–9566 (2017).',
    'stow0817': 'Stow, S. M. et al. An Interlaboratory Evaluation of Drift Tube Ion Mobility–Mass Spectrometry Collision Cross Section Measurements. Anal. Chem. 89, 9048–9055 (2017).',
    'zhen0917': 'Zheng, X. et al. A structural examination and collision cross section database for over 500 metabolites and xenobiotics using drift tube ion mobility spectrometry. Chem. Sci. 8, 7724–7736 (2017).',
    'righ0218': 'Righetti, L. et al. Ion mobility-derived collision cross section database: Application to mycotoxin analysis. Analytica Chimica Acta 1014, 50–57 (2018).',
    'pagl0314': 'Paglia, G. et al. Ion Mobility Derived Collision Cross Sections to Support Metabolomics Applications. Anal. Chem. 86, 3985–3993 (2014).',
    'moll0218': 'Mollerup, C. B., Mardal, M., Dalsgaard, P. W., Linnet, K. & Barron, L. P. Prediction of collision cross section and retention time for broad scope screening in gradient reversed-phase liquid chromatography-ion mobility-high resolution accurate mass spectrometry. Journal of Chromatography A 1542, 82–88 (2018).',
    'nich1118': 'Nichols, C. M. et al. Untargeted Molecular Discovery in Primary Metabolism: Collision Cross Section as a Molecular Descriptor in Ion Mobility-Mass Spectrometry. Anal. Chem. 90, 14484–14492 (2018).',
    'may_0114': 'May, J. C. et al. Conformational Ordering of Biomolecules in the Gas Phase: Nitrogen Collision Cross Sections Measured on a Prototype High Resolution Drift Tube Ion Mobility-Mass Spectrometer. Anal. Chem. 86, 2107–2116 (2014).',
    'hine0817': 'Hines, K. M., Ross, D. H., Davidson, K. L., Bush, M. F. & Xu, L. Large-Scale Structural Characterization of Drug and Drug-Like Compounds by High-Throughput Ion Mobility-Mass Spectrometry. Anal. Chem. 89, 9023–9030 (2017).',
    'hine1217': 'Hines, K. M. et al. Characterization of the Mechanisms of Daptomycin Resistance among Gram-Positive Bacterial Pathogens by Multidimensional Lipidomics. mSphere 2, 99–16 (2017).',
    'hine0217': 'Hines, K. M., Herron, J. & Xu, L. Assessment of altered lipid homeostasis by HILIC-ion mobility-mass spectrometry-based lipidomics. The Journal of Lipid Research 58, 809–819 (2017).',
    'groe0815': 'Groessl, M., Graf, S. & Knochenmuss, R. High resolution ion mobility-mass spectrometry for separation and identification of isomeric lipids. Analyst 140, 6904–6911 (2015).',
    'bijl0517': 'Bijlsma, L. et al. Prediction of Collision Cross-Section Values for Small Molecules: Application to Pesticide Residue Analysis. Anal. Chem. 89, 6583–6589 (2017).',
    'hine0119': 'Hines, K. M. & Xu, L. Lipidomic consequences of phospholipid synthesis defects in Escherichia coli revealed by HILIC-ion mobility-mass spectrometry. Chemistry and Physics of Lipids 219, 15–22 (2019).',
    'leap0219': 'Leaptrot, K. L., May, J. C., Dodds, J. N. & McLean, J. A. Ion mobility conformational lipid atlas for high confidence lipidomics. Nature Communications 1–9 (2019).',
    'blaz0818': 'Blaženović, I. et al. Increasing Compound Identification Rates in Untargeted Lipidomics Research with Liquid Chromatography Drift Time–Ion Mobility Mass Spectrometry. Anal. Chem. 90, 10758–10764 (2018).',
    'vasi0120': 'Vasilopoulou, C. G. et al. Trapped ion mobility spectrometry and PASEF enable in-depth lipidomics from minimal sample amounts. Nature Communications 1–11 (2020).',
    'tsug0220': 'Tsugawa, H. et al. MS-DIAL 4: accelerating lipidomics using an MS/MS, CCS, and retention time atlas. bioRxiv 37, 513 (2020).',
    'lian0118': 'Lian, R. et al. Ion mobility derived collision cross section as an additional measure to support the rapid analysis of abused drugs and toxic compounds using electrospray ion mobility time-of-flight mass spectrometry. Anal. Methods 10, 749–756 (2018).',
    'teja0918': 'Tejada-Casado, C. et al. Collision cross section (CCS) as a complementary parameter to characterize human and veterinary drugs. Analytica Chimica Acta 1043, 52–63 (2018).',
    'pola0620': 'Poland, J. C. et al. Collision Cross Section Conformational Analyses of Bile Acids via Ion Mobility–Mass Spectrometry. Journal of the American Society for Mass Spectrometry 31, 1625–1631 (2020).',
    'dodd0220': 'Dodds, J. et al. Rapid Characterization of Per- and Polyfluoroalkyl Substances (PFAS) by Ion Mobility Spectrometry−Mass Spectrometry (IMS-MS). Anal. Chem. 92, 4427-4435 (2020).',
    'celm1120': 'Celma, A. et al. Improving Target and Suspect Screening High-Resolution Mass Spectrometry Workflows in Environmental Analysis by Ion Mobility Separation. Environ. Sci. Technol. 54, 15120-15131 (2020)'
}


# ordered lists of src_tags and thier associated full reference info
# the index of each corresponds to the reference number - 1 
#   i.e. src_tags_ordered[0] -> src_tag for reference 1
#   and likewise references_ordered[0] -> full reference info for reference 1
src_tags_ordered = [
    'may_0114',
    'pagl0314',
    'groe0815',
    'zhou1016',
    'hine0217',
    'bijl0517',
    'hine0817',
    'stow0817',
    'zhou0817',
    'zhen0917',
    'hine1217',
    'lian0118',
    'moll0218',
    'righ0218',
    'teja0918',
    'nich1118',
    'hine0119',
    'leap0219',
    'blaz0818',
    'vasi0120',
    'tsug0220',
    'pola0620',
    'dodd0220',
    'celm1120'
]

references_ordered = [
    'May, J. C. et al. Conformational Ordering of Biomolecules in the Gas Phase: Nitrogen Collision Cross Sections Measured on a Prototype High Resolution Drift Tube Ion Mobility-Mass Spectrometer. Anal. Chem. 86, 2107–2116 (2014).',
    'Paglia, G. et al. Ion Mobility Derived Collision Cross Sections to Support Metabolomics Applications. Anal. Chem. 86, 3985–3993 (2014).Paglia, G. et al. Ion Mobility Derived Collision Cross Sections to Support Metabolomics Applications. Anal. Chem. 86, 3985–3993 (2014).',
    'Groessl, M., Graf, S. & Knochenmuss, R. High resolution ion mobility-mass spectrometry for separation and identification of isomeric lipids. Analyst 140, 6904–6911 (2015).',
    'Zhou, Z., Shen, X., Tu, J. & Zhu, Z.-J. Large-Scale Prediction of Collision Cross-Section Values for Metabolites in Ion Mobility-Mass Spectrometry. Anal. Chem. 88, 11084–11091 (2016).',
    'Hines, K. M., Herron, J. & Xu, L. Assessment of altered lipid homeostasis by HILIC-ion mobility-mass spectrometry-based lipidomics. The Journal of Lipid Research 58, 809–819 (2017).',
    'Bijlsma, L. et al. Prediction of Collision Cross-Section Values for Small Molecules: Application to Pesticide Residue Analysis. Anal. Chem. 89, 6583–6589 (2017).',
    'Hines, K. M., Ross, D. H., Davidson, K. L., Bush, M. F. & Xu, L. Large-Scale Structural Characterization of Drug and Drug-Like Compounds by High-Throughput Ion Mobility-Mass Spectrometry. Anal. Chem. 89, 9023–9030 (2017).',
    'Stow, S. M. et al. An Interlaboratory Evaluation of Drift Tube Ion Mobility–Mass Spectrometry Collision Cross Section Measurements. Anal. Chem. 89, 9048–9055 (2017).',
    'Zhou, Z., Tu, J., Xiong, X., Shen, X. & Zhu, Z.-J. LipidCCS: Prediction of Collision Cross-Section Values for Lipids with High Precision To Support Ion Mobility–Mass Spectrometry-Based Lipidomics. Anal. Chem. 89, 9559–9566 (2017).',
    'Zheng, X. et al. A structural examination and collision cross section database for over 500 metabolites and xenobiotics using drift tube ion mobility spectrometry. Chem. Sci. 8, 7724–7736 (2017).',
    'Hines, K. M. et al. Characterization of the Mechanisms of Daptomycin Resistance among Gram-Positive Bacterial Pathogens by Multidimensional Lipidomics. mSphere 2, 99–16 (2017).',
    'Lian, R. et al. Ion mobility derived collision cross section as an additional measure to support the rapid analysis of abused drugs and toxic compounds using electrospray ion mobility time-of-flight mass spectrometry. Anal. Methods 10, 749–756 (2018).',
    'Mollerup, C. B., Mardal, M., Dalsgaard, P. W., Linnet, K. & Barron, L. P. Prediction of collision cross section and retention time for broad scope screening in gradient reversed-phase liquid chromatography-ion mobility-high resolution accurate mass spectrometry. Journal of Chromatography A 1542, 82–88 (2018).',
    'Righetti, L. et al. Ion mobility-derived collision cross section database: Application to mycotoxin analysis. Analytica Chimica Acta 1014, 50–57 (2018).',
    'Tejada-Casado, C. et al. Collision cross section (CCS) as a complementary parameter to characterize human and veterinary drugs. Analytica Chimica Acta 1043, 52–63 (2018).',
    'Nichols, C. M. et al. Untargeted Molecular Discovery in Primary Metabolism: Collision Cross Section as a Molecular Descriptor in Ion Mobility-Mass Spectrometry. Anal. Chem. 90, 14484–14492 (2018).',
    'Hines, K. M. & Xu, L. Lipidomic consequences of phospholipid synthesis defects in Escherichia coli revealed by HILIC-ion mobility-mass spectrometry. Chemistry and Physics of Lipids 219, 15–22 (2019).',
    'Leaptrot, K. L., May, J. C., Dodds, J. N. & McLean, J. A. Ion mobility conformational lipid atlas for high confidence lipidomics. Nature Communications 1–9 (2019).',
    'Blaženović, I. et al. Increasing Compound Identification Rates in Untargeted Lipidomics Research with Liquid Chromatography Drift Time–Ion Mobility Mass Spectrometry. Anal. Chem. 90, 10758–10764 (2018).'
    'Vasilopoulou, C. G. et al. Trapped ion mobility spectrometry and PASEF enable in-depth lipidomics from minimal sample amounts. Nature Communications 1–11 (2020).',
    'Tsugawa, H. et al. MS-DIAL 4: accelerating lipidomics using an MS/MS, CCS, and retention time atlas. bioRxiv 37, 513 (2020).',
    'Poland, J. C. et al. Collision Cross Section Conformational Analyses of Bile Acids via Ion Mobility–Mass Spectrometry. Journal of the American Society for Mass Spectrometry 31, 1625–1631 (2020).',
    'Dodds, J. et al. Rapid Characterization of Per- and Polyfluoroalkyl Substances (PFAS) by Ion Mobility Spectrometry−Mass Spectrometry (IMS-MS). Anal. Chem. 92, 4427-4435 (2020).',
    'Celma, A. et al. Improving Target and Suspect Screening High-Resolution Mass Spectrometry Workflows in Environmental Analysis by Ion Mobility Separation. Environ. Sci. Technol. 54, 15120-15131 (2020)'
]

reference_links_ordered = [
    'https://pubs.acs.org/doi/abs/10.1021/ac4038448',
    'https://pubs.acs.org/doi/abs/10.1021/ac500405x',
    'https://pubs.rsc.org/en/content/articlelanding/2015/an/c5an00838g',
    'https://pubs.acs.org/doi/abs/10.1021/acs.analchem.6b03091',
    'http://www.jlr.org/content/58/4/809.abstract',
    'https://pubs.acs.org/doi/abs/10.1021/acs.analchem.7b00741',
    'https://pubs.acs.org/doi/abs/10.1021/acs.analchem.7b01709',
    'https://pubs.acs.org/doi/abs/10.1021/acs.analchem.7b01729',
    'https://pubs.acs.org/doi/abs/10.1021/acs.analchem.7b02625',
    'https://pubs.rsc.org/en/content/articlelanding/2017/sc/c7sc03464d',
    'https://msphere.asm.org/content/2/6/e00492-17',
    'https://pubs.rsc.org/en/content/articlelanding/2018/ay/c7ay02808c',
    'https://www.sciencedirect.com/science/article/pii/S0021967318301894',
    'https://www.sciencedirect.com/science/article/pii/S0003267018301521',
    'https://www.sciencedirect.com/science/article/abs/pii/S0003267018311693',
    'https://pubs.acs.org/doi/abs/10.1021/acs.analchem.8b04322',
    'https://www.sciencedirect.com/science/article/pii/S0009308418302238',
    'https://www.nature.com/articles/s41467-019-08897-5',
    'https://pubs.acs.org/doi/10.1021/acs.analchem.8b01527',
    'https://www.nature.com/articles/s41467-019-14044-x',
    'https://www.biorxiv.org/content/10.1101/2020.02.11.944900v1',
    'https://pubs.acs.org/doi/10.1021/jasms.0c00015',
    'https://pubs.acs.org/doi/10.1021/acs.analchem.9b05364',
    'https://pubs.acs.org/doi/10.1021/acs.est.0c05713'
]



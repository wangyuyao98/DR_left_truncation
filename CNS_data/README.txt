This repository contains the CNS lymphoma data ("CNS_data.txt") and the R code ("CNS_analysis.R") used for the data analyses in the paper. 
Permission is granted to anyone wishing to use this data set which was derived from a study of methotrexate-based chemotherapy (Wang et al., 2015).

Variable description:
relapse = 1 if relapse was observed, 0 otherwise
PFS = time-to-1st-relapse if relapse=1, otherwise time-to-end-of-follow-up
relapse2 = 1 if second relapse was observed, 0 otherwise
PFS2 = time-to-2nd-relapse if relapse2=1, otherwise time-to-end-of-follow-up
death = 1 if death was observed, 0 otherwise
OS = time-to-death if death=1, otherwise time-to-end-of-follow-up
gender: 0=women, 1=men
KPS_diag = Karnofsky performance score at diagnosis 
site = Initial site: 1=brain only, 2=cerebrospinal fluid(csf), 3=eye
chemo = Initial chemotherapy: 1=methotrexate(MTX) only, 2=MTX polychemo, 3=other
radiation = Initial radiation therapy: 1=Yes, 0=No

References:
Wang, N., Gill, C., Betensky, R. A., and Batchelor, T. (April, 2015). Relapse patterns in primary CNS diuse large B-cell lymphoma. Annual Meeting of the American Academy of Neurology, Washington D.C.

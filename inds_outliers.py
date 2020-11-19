#rat1Ind = {'230': 3, '355': 2, '358': 2, '362': 2, '673': 2, '674': 2, '816': 3, '825': 3}
rat1Ind = {'230': 3, '355': 2, '358': None, '362': 2, '673': 2, '674': 2, '816': 3, '825': 3}
shockInd = {'230': 5, '355': 5, '358': 5, '362': 5, '673': 5, '674': 5, '816': 5, '825': 5}
rat2Ind = {'230': 9, '355': 3, '358': 3, '362': 3, '673': 3, '674': 3, '816': None, '825': None}
epmInd = {'230': 1, '355': 1, '358': 1, '362': 1, '673': 1, '674': 1, '816': 1, '825': 1}
shockExtInd = {'230': 6, '355': 6, '358': 6, '362': 6, '673': None, '674': None, '816': 6, '825': 6}
toyRatInd = {'230': 2, '355': 7, '358': 7, '362': None, '673': 6, '674': 6, '816': 2, '825': 2}
#shockHabInd = {'230': 4, '355': 8, '358': 8, '362': 8, '673': 7, '674': 7, '816': 7, '825': None}
shockHabInd = {'230': 4, '355': 8, '358': 8, '362': 8, '673': 7, '674': 7, '816': 7, '825': 4}

outliersEPM = {'230': 3, '355': 2, '358': 3, '362': 1, '673': 3, '674': 2, '816': 7, '825': 4}
outliersRat1 = {'230': 5, '355': 2, '358': 2, '362': 1, '673': 0, '674': 0, '816': 5, '825': 0}
outliersShock = {'230': 3, '355': 2, '358': 2, '362': 1, '673': 1, '674': 1, '816': 1, '825': 5}
outliersRat2 = {'230': 4, '355': 0, '358': 1, '362': 1, '673': 2, '674': 1, '816': None, '825': None}
outliersShockExt = {'230': 1, '355': 0, '358': 2, '362': 2, '673': None, '674': None, '816': 1, '825': 5}
outliersToyRat = {'230': 4, '355': 0, '358': 4, '362': 0, '673': 0, '674': 0, '816': 6, '825': 6}
outliersShockHab =  {'230': 0, '355': 0, '358': 2, '362': 0, '673': 1, '674': 0, '816': 2, '825': 0}

assays = ['epm', 'rat1', 'rat2', 'shock', 'shockext', 'toyrat', 'shockhab']
mice = ['230', '355', '358', '362', '673', '674', '816', '825']

Inds = {}
Inds['rat1'] = rat1Ind
Inds['shock'] = shockInd
Inds['rat2'] = rat2Ind
Inds['epm'] = epmInd
Inds['shockext'] = shockExtInd
Inds['toyrat'] = toyRatInd
Inds['shockhab'] = shockHabInd

outliers = {}
outliers['epm'] = outliersEPM
outliers['rat1'] = outliersRat1
outliers['shock'] = outliersShock
outliers['rat2'] = outliersRat2
outliers['shockext'] = outliersShockExt
outliers['toyrat'] = outliersToyRat
outliers['shockhab'] = outliersShockHab

## denote the pcas for which there is weird noise
rmpcas = {}
for assay in assays:
    rmpcas[assay] = {}
    for mouse in mice:
        rmpcas[assay][mouse] = []
rmpcas['rat1']['230'] = [2]
rmpcas['rat1']['358'] = [1]
rmpcas['shock']['362'] = [1, 2]
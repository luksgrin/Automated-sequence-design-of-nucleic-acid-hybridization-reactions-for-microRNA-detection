import sys
import subprocess
import random
import datetime
import time
#ViennaRNA python3 library (https://www.tbi.univie.ac.at/RNA/documentation.html)
sys.path.append('/home/lugoibel/ViennaRNA/interfaces/Python3')
#NUPACK suite (http://www.nupack.org/partition/new) that employs a modified wrapper originally by:
#Salis, H., Mirsky, E., & Voigt, C. (2009).
#Automated design of synthetic ribosome binding sites to control protein expression.
#Nature Biotechnology, 27(10), 946-950.
sys.path.append('/home/lugoibel/nupack3.2.2/python')
import RNA
from NuPACK import NuPACK
import plotly.plotly as py
import plotly.offline as offline
import plotly.graph_objs as go

#########################################
#                                       #
#       DEFINITION OF PARAMETERS        #
#                                       #
#########################################

#Start time
start_time = time.time()
#Vienna parameters:
    #Mathews parameterfile
RNA.read_parameter_file(
    '/home/lugoibel/ViennaRNA/misc/dna_mathews2004.par')
    #No dangles
RNA.cvar.dangles = 0
    #No coversion from DNA into RNA
RNA.cvar.nc_fact = 1
#Define nucleotides
NUCS = ['A','T','G','C']
#Define circuit sequence names
GUIDE = [
    'miRNA',
    'sensor',
    'transducer',
    'clamp',
    'T7p',
    'fuel']
#Boltzmann function parameters
BETA = 1/0.593
NUM_e = 2.7182818284590452353
DGbp = -1.25
#Metropolis parameters
Bm0 = 1100        #con 1e3  no converge
D = 1.00007
#Define shadow circuit constant components
shdw = {'S2': 'TGAGATGTAAAGGATGAGTGAGATG',
    'T2': 'CACTCATCCTTTACATCTCAAACACTCTATTCA'}

#########################################
#                                       #
#       DEFINITION OF FUNCTIONS         #
#                                       #
#########################################

#Define command line input system
def cmdinput():
    global USERINPUT
    global GUIDE
    looping = True
    while looping:
        if 'U' in USERINPUT:
            USERINPUT = USERINPUT.replace(
                'U','T')
        UNIQ = set(USERINPUT)
        #Checks if input is a sequence of adequate length
        if (UNIQ.issubset(NUCS) and
                len(USERINPUT) >= 20):
            seqs_preit[GUIDE[0]] = USERINPUT[:25]
            looping = False
        #Checks if input is meant to be a test
        elif USERINPUT == 'TEST':
            GUIDE = ['Rodrigo_miRNA'] + GUIDE[1:]
            seqs_preit['Rodrigo_miRNA'] = 'TGGAGTGTGACAATGGTGTTTG'
            looping = False
        #Exit system
        elif USERINPUT == 'EXIT':
            exit()
        #Retry input if previous statements are false
        else:
            USERINPUT = input(
                'Enter a VALID input: '
                ).upper()
    return None

#Define fasta file input system. Saves data in a dictionary as
#key = header and value = sequence, only if the sequence is
#adequate
def fileinput():
    dict = {}

    for line in open(USERINPUT):
        line = line.strip('\n')

        if line[0] == '>':
            key = line[1:].split()[0]
            value = ''

        else:
            value += line

        if (set(value).issubset(NUCS) and
                len(value) >= 20):
            dict[key] = value[:25]
    return dict

#Define reverse complementary generator
def revcomp(seq):
    seq = seq.upper(
        ).replace('A','t'
        ).replace('T','a'
        ).replace('G','c'
        ).replace('C','g'
        )[::-1].upper()
    return seq

#Random sequence builder
def randseq(length):
    out = ''
    for n in range(length):
        out += random.sample(NUCS, 1)[0]
    return out

#Define circuit core sequences generator
def genseq(miRNA, prom):
    n = len(miRNA)
    rootseq = (miRNA.upper()
        + randseq(5)#'TATTC'
        + prom)
    sensor = revcomp(
        rootseq[: n+8])
    transducer = rootseq[6:]
    clamp = revcomp(
        rootseq[n - 1 :])
    fuel = rootseq[6: n + 8]
    return (sensor,
        transducer,
        clamp,
        fuel)

#Define Boltzmann function
def bolfunc(seq1, seq2, seq_DG):
    Pairkey = (seq1
        + '_'
        + seq2)
    
    Numerator = NUM_e**(- BETA*seq_DG[Pairkey])
    Denominator = Numerator

    if seq1 == GUIDE[0]:
        SecondKey = 'sensor_transducer'

    elif seq1 == 'transducer':
        SecondKey = 'clamp_T7p'

    elif seq1 == 'fuel':
        SecondKey = GUIDE[0] + '_sensor'

    Denominator += NUM_e**(- BETA*seq_DG[SecondKey])
    func = Numerator/Denominator

    return func

#Define function for probability calculation employed in
#secondary pairments
def probfunc(seq1, seq2, seq_DG, seqs):
    Pairkey = (seq1
        + '_'
        + seq2)
    Numerator = NUM_e**(- BETA*seq_DG[Pairkey])

    if seq1 == 'sensor':
        L = 19

    if seq1 == 'clamp':
        L = len(seqs[GUIDE[4]])
    Denominator = NUM_e**(- BETA*L*DGbp)
    func = Numerator/Denominator

    if func > 1:
        func = 1
    return func

#Define toehold score function
def toeholdscore(name, seq_ss):
    DIST = (len(seqs_preit['transducer'])
        - len(seqs_preit['T7p'])
        + 3)
    struct = seq_ss[name].split(
        '&'
        )[1][(DIST-6):DIST]
    j = 0

    for symbol in struct:
        if symbol == '.':
            j += 1
    return j

#Define Packing and Scoring function.
def scorefunc(seqs):
    seq_DG = {}
    seq_ss = {}
    i = -1
    #Saves in a dictionary the MFE and structure of circuit pairs
    for seq1 in GUIDE[:-2]:
        i += 1
        seq2 = GUIDE[i + 1]
        name = (seq1
            + '_'
            + seq2)
        (ss, mfe) = RNA.cofold(
            (seqs[seq1]
            + '&'
            + seqs[seq2]))
        seq_DG[name] = mfe
        seq_ss[name] = (ss[: len(seqs[seq1])]
            + '&'
            + ss[(len(seqs[seq1])) :-1])

    (ss, mfe) = RNA.cofold(
        (seqs['fuel']
        + '&'
        + seqs['sensor']))
    seq_DG['fuel_sensor'] = mfe
    seq_ss['fuel_sensor'] = (ss[: len(seqs['fuel'])]
        + '&'
        + ss[(len(seqs['fuel'])) :-1])
    #Caulculates pair probabilities and Score

    P1 = bolfunc(
        GUIDE[0],
        'sensor',
        seq_DG)
    P2 = bolfunc(
        'transducer',
        'clamp',
        seq_DG)
    P3 = probfunc(
        'sensor',
        'transducer',
        seq_DG,
        seqs)
    P4 = probfunc(
        'clamp',
        'T7p',
        seq_DG,
        seqs)
    P5 = bolfunc(
        'fuel',
        'sensor',
        seq_DG)
    T = toeholdscore('sensor_transducer', seq_ss)

    score = P1*P2*P3*P4*P5*(6-T)/6
    dats = [P1,P2,P3,P4,P5,T,score]
    return dats

#Define mutation function
def mutf(seqs):
    seqs_aftermutation = {}
    
    #Creates a new dictionary with sequences
    for element in seqs:
        seqs_aftermutation[element] = seqs[element]
    
    #Creates a new guidelist excluding miRNA and T7p
    mutlist = GUIDE[1:-2] + [GUIDE[-1]]
    
    #Chooses a random base from a random sequence from ensemble
    target_name = random.sample(mutlist, 1)[0]
    target_seq = list(seqs[target_name])
    position = random.randint(0, (len(target_seq) - 1))
    base = random.sample(NUCS, 1)[0]

    while base == target_seq[position]:
        base = random.sample(NUCS, 1)[0]
    
    #Writes the mutated sequence
    target_seq[position] = base
    target_seq = ''.join(target_seq)
    seqs_aftermutation[target_name] = target_seq
    
    return seqs_aftermutation

#Define a function that interprets NuPACK output files
def eqcon(dict, guide):
    outlist = []
    outdict = {}

    for el in dict['complexes_concentrations']:
        stand = round((float(el[-1])/1e-8), 2)

        if stand < 0.1:
            continue

        cmplx = list(map(int, el[0:-2]))

        name = []
        i = -1
        for n in cmplx:
            i += 1

            if n:
                name += n*[guide[i]]

        name = '_'.join(name)

        outlist += [name]
        outdict[name] = [el[-1], stand]
    return outlist, outdict

#Define test-tube prediction of final equilibriums by means of NuPACK
def test_tube(seqs, guide):
    print('Calculating test-tube NuPACK simulation')
    seq_list = []
    concent = [1e-6, 1e-6]

    if 'fuel' not in guide:
        concent += [1e-6, 1e-6]

    for el in guide:
        seq_list += [seqs[el]]

    eq_1 = NuPACK(
        Sequence_List=seq_list,
        material='dna')
    eq_2 = NuPACK(
        Sequence_List=seq_list,
        material='dna')

    eq_1.complexes(
        dangles='none',
        MaxStrands=2,
        quiet=True)
    eq_2.complexes(
        dangles='none',
        MaxStrands=2,
        quiet=True)

    eq_1.concentrations(
        concentrations=[1e-6] + concent,
        quiet=True)
    eq_2.concentrations(
        concentrations=[1e-9] + concent,
        quiet=True)

    (eq_1order, eq_1) = eqcon(eq_1, guide)
    (eq_2order, eq_2) = eqcon(eq_2, guide)
    EQUILIBRIUMGUIDES = [eq_1order, eq_2order]
    return EQUILIBRIUMGUIDES, eq_1, eq_2

#Define bar-chart plot function for NuPACK test-tube prediction
def eqsbarplot(guides, dict1, dict2):
    global timessufix
    dat1 = []
    dat2 = []

    for list in guides:
        for el in list:

            if guides[0] == list:
                dat1 += [dict1[el][-1]]

            else:
                dat2 += [dict2[el][-1]]

    trace1 = go.Bar(
        x=guides[0],
        y=dat1,
        name='With input')
    trace2 = go.Bar(
        x=guides[1],
        y=dat2,
        name='Without input')

    data = [trace1, trace2]
    layout = go.Layout(
        barmode='group',
        title='Equilibrium concentrations for species',
        yaxis=dict(title='% abundance'))
    fig = go.Figure(
        data=data,
        layout=layout)
    filename = ('Equilibrium_study_'
        + timesuffix
        + '.html')
    offline.plot(
        fig,
        filename=filename,
        auto_open=False)

    return None

#Define metropolis function to induce random sampling
def Metropolis():
    global Dats_preit, Score_preit, seqs_preit
    Bmk = Bm0*(D**k)
    M = NUM_e**(
        - Bmk*(
            Score_preit
            - Score_posit))

    if random.random() < M:
#       print('\nMetropolis MUTATED\n')
        Dats_preit = Dats_posit
        Score_preit = Score_posit
        seqs_preit = seqs_posit
    return None

#Define percentage progress percentage function
def progress():
    global perc_0
    perc_1 = (k/100000)*100

    if int(perc_1/5) > int(perc_0/5):
        perc_0 = perc_1
        print(
            'Status: '
            + str(int(perc_0))
            + '% completed')
    return None

#Define shadow circuit generation function
def shadowcirc(transducer):
    outdict = {}
    
    for el in shdw:
        outdict[el] = shdw[el]
    
    MFE = RNA.cofold(
        seqs_preit['sensor']
        + '&'
        + seqs_preit['transducer'])[1]

    mfe = RNA.cofold(
        outdict['S2']
        + '&'
        + outdict['T2'])[1]
    
    b_area = outdict['S2'][:-5]
    
    if MFE < -31:
        times = int((MFE + 31)/3) + 4
        
        for n in range(times):
            b_area += random.sample(NUCS, 1)[0]
    i = 0
    while abs(MFE - mfe) > 0:
        
        i += 1
        target_index = random.randint(0, (len(b_area) - 1))

        b_area = list(b_area)
        base = random.sample(NUCS, 1)[0]

        while base == b_area[target_index]:
            base = random.sample(NUCS, 1)[0]

        b_area[target_index] = base
        b_area = ''.join(b_area)

        outdict['S2'] = (b_area
            + outdict['S2'][-5:])

        outdict['T2'] = (revcomp(b_area)
            + outdict['S2'][-13:])

        mfe = RNA.cofold(
            outdict['S2']
            + '&'
            + outdict['T2'])[1]

        if i == 1000:
            break
        
        for el in NUCS:
            if (4*el) in b_area:
                mfe = 1e3
        
    master = (outdict['T2'][-20:]
        + transducer[:19])
    
    AND_clamp = master[7:-6]
    AND = revcomp(master)
    
    outdict['AND_clamp'] = AND_clamp
    outdict['AND'] = AND
    
    keyss = []
    for el in outdict.keys():
        keyss += [el]
    keyss.sort()
    
    return outdict, keyss 

#MAIN
def main():
    global k, timesuffix, perc_0, seqs_preit, seqs_posit
    global Score_preit, Score_posit, Dats_preit, Dats_posit
        
    #Moment in time:
    timesuffix = '_'.join(
        str(datetime.datetime.now()
        ).split())

    (seqs_preit['sensor'],
    seqs_preit['transducer'],
    seqs_preit['clamp'],
    seqs_preit['fuel']) = genseq(seqs_preit[GUIDE[0]], seqs_preit['T7p'])

    Dats_preit = scorefunc(seqs_preit)
    Score_preit = Dats_preit[-1]

    (equilibriumguide,
    eq_1,
    eq_2) = test_tube(seqs_preit, GUIDE[:-1])

    fuelguide = GUIDE[:2] + [GUIDE[-1]]
    fuelguide = fuelguide[::-1]

    (equilibriumguide_fuel,
    w_fuel,
    wo_fuel) = test_tube(seqs_preit, fuelguide)
        
    OUTFILE = open(
        'Output_'
        + GUIDE[0]
        + '_'
        + timesuffix
        + '.txt',
        'w')
    OUTFILE.write('This is the output of your job done on '
        + timesuffix
        + '\n')

    for el in GUIDE:
        OUTFILE.write('>'
            + el
            + '\n'
            + seqs_preit[el]
            + '\n')
        
    OUTFILE.write('\nP1 = ' + str(Dats_preit[0]) + '\n')
    OUTFILE.write('P2 = ' + str(Dats_preit[1]) + '\n')
    OUTFILE.write('P3 = ' + str(Dats_preit[2]) + '\n')
    OUTFILE.write('P4 = ' + str(Dats_preit[3]) + '\n')
    OUTFILE.write('P5 = ' + str(Dats_preit[4]) + '\n')
    OUTFILE.write('Toehold = ' + str(Dats_preit[5]) + '\n')
    OUTFILE.write('Score = ' + str(Score_preit) + '\n')
    OUTFILE.write('Standarized score = '
        + str(Score_preit*100/Dats_preit[3])
        + '\n')


    OUTFILE.write('\n------WITH INPUT------')
    OUTFILE.write('\nComplexes')
    OUTFILE.write('\tConcentration (M)')
    OUTFILE.write('\tStandarized (%)\n')

    for el in equilibriumguide[0]:
        OUTFILE.write(el
            + '\t'
            + eq_1[el][0]
            + '\t'
            + str(eq_1[el][1])
            + '\n')

    OUTFILE.write('\n------WITHOUT INPUT------')
    OUTFILE.write('\nComplexes')
    OUTFILE.write('\tConcentration (M)')
    OUTFILE.write('\tStandarized (%)\n')

    for el in equilibriumguide[1]:
        OUTFILE.write(el
            + '\t'
            + eq_2[el][0]
            + '\t'
            + str(eq_2[el][1])
            + '\n')

    OUTFILE.write('\nFuel transduction assessment\n')
    OUTFILE.write('\n------WITH FUEL------')
    OUTFILE.write('\nComplexes')
    OUTFILE.write('\tConcentration (M)')
    OUTFILE.write('\tStandarized (%)\n')

    for el in equilibriumguide_fuel[0]:
        OUTFILE.write(el
            + '\t'
            + w_fuel[el][0]
            + '\t'
            + str(w_fuel[el][1])
            + '\n')

    OUTFILE.write('\n------WITHOUT FUEL------')
    OUTFILE.write('\nComplexes')
    OUTFILE.write('\tConcentration (M)')
    OUTFILE.write('\tStandarized (%)\n')

    for el in equilibriumguide_fuel[1]:
        OUTFILE.write(el
            + '\t'
            + wo_fuel[el][0]
            + '\t'
            + str(wo_fuel[el][1])
            + '\n')
    OUTFILE.write('\n')

    #1e5 cycles of mutations and selection following the global score
    
    k = 0
    perc_0 = 0

    for n in range(int(1e5)):
        k += 1
        seqs_posit = mutf(seqs_preit)
        Dats_posit = scorefunc(seqs_posit)
        Score_posit = Dats_posit[-1]

        if Score_posit >= Score_preit:
            Dats_preit = Dats_posit
            Score_preit = Score_posit
            seqs_preit = seqs_posit

        else:
            Metropolis()
            
        progress()
        
    (equilibriumguide,
    eq_1,
    eq_2) = test_tube(seqs_preit, GUIDE[:-1])
    eqsbarplot(equilibriumguide,
        eq_1,
        eq_2)

    (equilibriumguide_fuel,
    w_fuel,
    wo_fuel) = test_tube(seqs_preit, fuelguide)
    #OUTFILE = open('Output_'+GUIDE[0]+timesuffix+'.txt', 'w')
    #OUTFILE.write('This is the output of your job done on '+timesuffix+'\n')

    for el in GUIDE:
        OUTFILE.write('>'
            + el
            + '\n'
            + seqs_preit[el]
            + '\n')

    OUTFILE.write('\nP1 = ' + str(Dats_preit[0]) + '\n')
    OUTFILE.write('P2 = ' + str(Dats_preit[1]) + '\n')
    OUTFILE.write('P3 = ' + str(Dats_preit[2]) + '\n')
    OUTFILE.write('P4 = ' + str(Dats_preit[3]) + '\n')
    OUTFILE.write('P5 = ' + str(Dats_preit[4]) + '\n')
    OUTFILE.write('Toehold = ' + str(Dats_preit[5]) + '\n')
    OUTFILE.write('Score = ' + str(Score_preit) + '\n')
    OUTFILE.write('Standarized score = '
        + str(Score_preit*100/Dats_preit[3])
        + '\n')

    OUTFILE.write('\n------WITH INPUT------')
    OUTFILE.write('\nComplexes')
    OUTFILE.write('\tConcentration (M)')
    OUTFILE.write('\tStandarized (%)\n')

    for el in equilibriumguide[0]:
        OUTFILE.write(el
            + '\t'
            + eq_1[el][0]
            + '\t'
            + str(eq_1[el][1])
            + '\n')

    OUTFILE.write('\n------WITHOUT INPUT------')
    OUTFILE.write('\nComplexes')
    OUTFILE.write('\tConcentration (M)')
    OUTFILE.write('\tStandarized (%)\n')

    for el in equilibriumguide[1]:
        OUTFILE.write(el
            + '\t'
            + eq_2[el][0]
            + '\t'
            + str(eq_2[el][1])
            + '\n')

    OUTFILE.write('\nFuel transduction assessment\n')
    OUTFILE.write('\n------WITH FUEL------')
    OUTFILE.write('\nComplexes')
    OUTFILE.write('\tConcentration (M)')
    OUTFILE.write('\tStandarized (%)\n')

    for el in equilibriumguide_fuel[0]:
        OUTFILE.write(el
            + '\t'
            + w_fuel[el][0]
            + '\t'
            + str(w_fuel[el][1])
            + '\n')

    OUTFILE.write('\n------WITHOUT FUEL------')
    OUTFILE.write('\nComplexes')
    OUTFILE.write('\tConcentration (M)')
    OUTFILE.write('\tStandarized (%)\n')

    for el in equilibriumguide_fuel[1]:
        OUTFILE.write(el
            + '\t'
            + wo_fuel[el][0]
            + '\t'
            + str(wo_fuel[el][1])
            + '\n')
    
    
    (shadow, shadowguide) = shadowcirc(
        seqs_preit['transducer'])
    
    OUTFILE.write('\nProposed shadow cancellation circuit\n')
    for el in shadowguide:
        OUTFILE.write('>'
            + el
            + '\n'
            + shadow[el]
            +'\n')

seqs_preit = {}

#T7p sequence
seqs_preit['T7p'] = 'GCGCTAATACGACTCACTATAGG'

#Define initial input
try:
    USERINPUT = sys.argv[1]
except:
    USERINPUT = input(
        'Enter your input: '
        ).upper()

#Checks if input is a raw sequence or a fasta file
if USERINPUT.lower().split('.')[-1] == 'fasta':
    insequences = fileinput()

    for el in insequences:
        GUIDE[0] = el
        seqs_preit[el] = insequences[el]
        main()

        for name in GUIDE[1:-1]:
            del seqs_preit[name]
else:
    cmdinput()
    main()

#NuPACK files cleanup
#subprocess.call(
#    'rm -r /home/lugoibel/nupack3.2.2/python/tmp*',
#    shell=True)

elapsed_time = str(
    (time.time() - start_time)/60)

print('Job finished on '
    + str(datetime.datetime.now()).split('.')[0]
    + '. Elapsed time was: '
    + elapsed_time[:-13]
    + ' minutes.')

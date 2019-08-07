import sys, json
from random import choice, random, randint
from datetime import datetime as date
from time import time
from math import exp
sys.path.append('/home/lugoibel/ViennaRNA/interfaces/Python3')
sys.path.append('/home/lugoibel/nupack3.2.2/python')
import RNA
cofold = RNA.cofold
from NuPACK import NuPACK, TmpCleaner
import plotly.plotly as py
import plotly.offline as offline
import plotly.graph_objs as go


#########################################
#                                       #
#       DEFINITION OF PARAMETERS        #
#                                       #
#########################################

#Start time
start_time = time()
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
#Boltzmann function parameters
BETA = 1/0.593
DGbp = -1.25
#Metropolis parameters
Bm0 = 1100        #con 1e3  no converge
D = 1.00007

#########################################
#                                       #
#        DEFINITION OF CLASSES          #
#                                       #
#########################################

#Define circuit class which contains circuit sequences and score
class Circuit:
    
    __GUIDE__ = ['miRNA', 'sensor', 'transducer', 'clamp', 'T7p', 'fuel']
    __mutlist__ = ['sensor', 'transducer', 'clamp', 'fuel']
    
    def __init__(self, name, miRNA):
        self.name = name
        self.miRNA = miRNA
        self.T7p = 'GCGCTAATACGACTCACTATAGG' #T7p sequence
        
        #Initial circuit sequences generator
        n = len(self.miRNA)
        rootseq = (self.miRNA.upper()
            + randseq(5)
            + self.T7p)
        self.sensor = revcomp(
            rootseq[: n+8])
        self.transducer = rootseq[6:]
        self.clamp = revcomp(
            rootseq[n - 1 :])
        self.fuel = rootseq[6: n+8]
        self.scoring()        
    
    #Duplicate class function
    def duplicate(self, circuit):
        for element in circuit.__dict__:
            setattr(self, element, getattr(circuit, element))
        
    #Scoring function
    def scoring(self):
        
        #Ideal Equilibrium complexes
        complexes = {'miRNA_sensor': Complex(self.miRNA, self.sensor),
            'sensor_transducer': Complex(self.sensor, self.transducer),
            'transducer_clamp': Complex(self.transducer, self.clamp),
            'clamp_T7p': Complex(self.clamp, self.T7p),
            'fuel_sensor': Complex(self.fuel, self.sensor)}
        
        #Define Boltzmann function
        def bolfunc(key):
            
            Numerator = exp(- BETA*complexes[key].DG)
            Denominator = Numerator

            if key == 'miRNA_sensor':
                SecondKey = 'sensor_transducer'

            elif key == 'transducer_clamp':
                SecondKey = 'clamp_T7p'

            elif key == 'fuel_sensor':
                SecondKey = 'miRNA_sensor'
            
            Denominator += exp(- BETA*complexes[SecondKey].DG)
            
            return Numerator/Denominator
        
        #Define function for probability calculation employed in
        #secondary pairments
        def probfunc(key):
            
            Numerator = exp(- BETA*complexes[key].DG)

            if key == 'sensor_transducer':
                L = 19

            if key == 'clamp_T7p':
                L = len(self.T7p)
            
            Denominator = exp(- BETA*L*DGbp)
            
            func = Numerator/Denominator

            if func > 1:
                func = 1
                
            return func
        
        #Define toehold score function
        def toeholdscore():
            
            DIST = (len(self.transducer)
                - len(self.T7p)
                + 3)
            struct = (complexes['sensor_transducer'].SS
                .split('&')
                [1]
                [(DIST-6):DIST])
            i = 0

            for symbol in struct:
                if symbol == '.':
                    i += 1
            
            return i
        
        #Calculates pair probabilities and Score
        self.score = (bolfunc('miRNA_sensor')
            *bolfunc('transducer_clamp')
            *probfunc('sensor_transducer')
            *probfunc('clamp_T7p')
            *bolfunc('fuel_sensor')
            *(6 - toeholdscore())/6)
        self.stdscore = self.score/probfunc('clamp_T7p')*100
        
        return None
    
    #Define mutation function
    def mutate(self):
                
        #Chooses a random base from a random sequence from ensemble
        target_name = choice(self.__mutlist__)
        target_seq = list(getattr(self, target_name))
        position = randint(0, (len(target_seq) - 1))
        base = choice(NUCS)

        while base == target_seq[position]:
            base = choice(NUCS)
    
        #Writes the mutated sequence
        target_seq[position] = base
        target_seq = ''.join(target_seq)
        setattr(self, target_name, target_seq)
        
        self.scoring()
        
        return None

#Define complex class
class Complex:
    
    def __init__(self, seq1, seq2):
        (ss, self.DG) = cofold(
            (seq1
            + '&'
            + seq2))
        self.SS = (ss[: len(seq1)]
            + '&'
            + ss[len(seq1) :-1]) 

#Define test-tube class with final equilibriums prediction by means of NuPACK
class Test_tube:

    def __init__(self, OBJECT, MODE='TOTAL'):
        self.name = OBJECT.name
        concent = [1e-6, 1e-6]
        
        if __name__ == '__main__':
            print('Calculating test-tube NuPACK simulation')
        
        if isinstance(OBJECT, Circuit) and MODE == 'TOTAL':
            guide = OBJECT.__GUIDE__[:-1]
            concent += [1e-6, 1e-6]
        elif isinstance(OBJECT, Circuit) and MODE == 'FUEL':
            guide = (OBJECT.__GUIDE__[0:2]
                + [OBJECT.__GUIDE__[-1]])[::-1]
        
        seq_list = []
        
        for el in guide:
            seq_list += [getattr(OBJECT, el)]
        
        eq1 = NuPACK(
            Sequence_List=seq_list,
            material='dna')
        eq2 = NuPACK(
            Sequence_List=seq_list,
            material='dna')
        
        eq1.complexes(
            dangles='none',
            MaxStrands=2,
            quiet=True)
        eq2.complexes(
            dangles='none',
            MaxStrands=2,
            quiet=True)
        
        eq1.concentrations(
            concentrations=[1e-6] + concent,
            quiet=True)
        eq2.concentrations(
            concentrations=[1e-9] + concent,
            quiet=True)
        
        (g1, eq1) = self.eqcon(eq1, guide)
        (g2, eq2) = self.eqcon(eq2, guide)
        
        self.eq1 = eq1
        self.eq2 = eq2
        self.__GUIDES__ = [g1, g2]
    
    #Define a function that interprets NuPACK output files
    def eqcon(self, dict, guide):
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


#Define shadow circuit generation function
class shadowcirc:
    
    __GUIDE__ = ['S2', 'T2', 'AND', 'AND_clamp']
    __mutlist__ = ['AND', 'AND_clamp']
        
    def __init__(self, circuitobject):
        self.S2 = 'TGAGATGTAAAGGATGAGTGAGATG'
        self.T2 = 'CACTCATCCTTTACATCTCAAACACTCTATTCA'    
            
        
        MFE = cofold(
            circuitobject.sensor
            + '&'
            + circuitobject.transducer)[1]

        mfe = cofold(
            self.S2
            + '&'
            + self.T2)[1]
    
        b_area = self.S2[:-5]
    
        if MFE < -31:
            times = int((MFE + 31)/3) + 4
        
            for n in range(times):
                b_area += choice(NUCS)
        i = 0
        while abs(MFE - mfe) > 0:
        
            i += 1
            target_index = randint(0, (len(b_area) - 1))

            b_area = list(b_area)
            base = choice(NUCS)

            while base == b_area[target_index]:
                base = choice(NUCS)

            b_area[target_index] = base
            b_area = ''.join(b_area)

            self.S2 = (b_area
            + self.S2[-5:])

            self.T2 = (revcomp(b_area)
                + self.S2[-13:])

            mfe = cofold(
                self.S2
                + '&'
                + self.T2)[1]

            if i == 1000:
                break
        
            for el in NUCS:
                if (4*el) in b_area:
                    mfe = 1e3
        
        master = (self.T2[-20:]
            + circuitobject.transducer[:19])
    
        self.AND_clamp = master[7:-6]
        self.AND = revcomp(master)
            
#########################################
#                                       #
#       DEFINITION OF FUNCTIONS         #
#                                       #
#########################################

#Define command line input system
def cmdinput():
    global USERINPUT, NAME, MIRNA
    looping = True
    while looping:
        if 'U' in USERINPUT:
            USERINPUT = USERINPUT.replace(
                'U','T')
        UNIQ = set(USERINPUT)
        #Checks if input is a sequence of adequate length
        if (UNIQ.issubset(NUCS) and
                len(USERINPUT) >= 20):
            NAME = 'miRNA_Circuit'
            MIRNA = USERINPUT[:25]
            looping = False
        #Checks if input is meant to be a test
        elif USERINPUT == 'TEST':
            NAME = 'Rodrigo_Circuit'
            MIRNA = 'TGGAGTGTGACAATGGTGTTTG'
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
            key = (line[1:].split()[0]
                + '_Circuit')
            value = ''

        else:
            value += line

        if (set(value).issubset(NUCS) and
                len(value) >= 20):
            dict[key] = value[:25]
    return dict

#Define reverse complementary generator
def revcomp(seq):
    seq = (seq
        .upper()
        .replace('A','t')
        .replace('T','a')
        .replace('G','c')
        .replace('C','g')
        [::-1]
        .upper())
    return seq

#Random sequence builder
def randseq(length):
    out = ''
    for n in range(length):
        out += choice(NUCS)
    return out

#Define bar-chart plot function for NuPACK test-tube prediction
def eqsbarplot(TEST_TUBE):
    global timessufix
    dat1 = []
    dat2 = []
    
    i = -1
    for list in TEST_TUBE.__GUIDES__:
        i +=1
        for el in list:
            if not i:
                dat1 += [TEST_TUBE.eq1[el][-1]]

            else:
                dat2 += [TEST_TUBE.eq2[el][-1]]

    trace1 = go.Bar(
        x=TEST_TUBE.__GUIDES__[0],
        y=dat1,
        name='With input')
    trace2 = go.Bar(
        x=TEST_TUBE.__GUIDES__[1],
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
    filename = ('Equilibrium_study_{}.html'
                    .format(timesuffix))
    offline.plot(
        fig,
        filename=filename,
        auto_open=False)

    return None

#Define metropolis function to induce random sampling
def Metropolis(ENSEMBLEpre, ENSEMBLEpost):
    Bmk = Bm0*(D**k)
    M = exp(
        - Bmk*(
            ENSEMBLEpre.score
            - ENSEMBLEpost.score))

    if random() < M:
        ENSEMBLEpre.duplicate(ENSEMBLEpost)
    
    return None

#Define data and header writing function
def OutfileWriter(OUTFILE,
                  ENSEMBLEpre,
                  SIMULATION1,
                  FUELSIMULATION1):

    def DataExtractor(TestTube, I):
        eqs = ['eq1', 'eq2']        
        for el in TestTube.__GUIDES__[I]:
            
            OUTFILE.write('{}\t{}\t{}\n'
                .format(el,
                        getattr(TestTube, eqs[I])[el][0],
                        getattr(TestTube, eqs[I])[el][1]))
            
        return None

    for el in ENSEMBLEpre.__GUIDE__:
        OUTFILE.write('>{}\n{}\n'
            .format(el, getattr(ENSEMBLEpre, el)))

    OUTFILE.write('Score = {}\nStandarized score = {}\n'
        .format(ENSEMBLEpre.score, ENSEMBLEpre.stdscore))
            
    OUTFILE.write('\n{0}WITH INPUT{0}'
        '\nComplexes\tConcentration (M)\tStandarized (%)\n'
        .format(6*'-'))

    DataExtractor(SIMULATION1, 0)
                  
    OUTFILE.write('\n{0}WITHOUT INPUT{0}'
        '\nComplexes\tConcentration (M)\tStandarized (%)\n'
        .format(6*'-'))
        
    DataExtractor(SIMULATION1, 1)
        
    OUTFILE.write('\nFuel transduction assessment\n'
        '\n{0}WITH FUEL{0}'
        '\nComplexes\tConcentration (M)\tStandarized (%)\n'
        .format(6*'-'))
    
    DataExtractor(FUELSIMULATION1, 0)

    OUTFILE.write('\n{0}WITHOUT FUEL{0}'
        '\nComplexes\tConcentration (M)\tStandarized (%)\n'
        .format(6*'-'))
        
    DataExtractor(FUELSIMULATION1, 1)
    
    OUTFILE.write('\n')
    
    return None

#Circuit object serialization function
def serialize_circuit(obj):    
    if isinstance(obj, Circuit):
        Ser = obj.__dict__
        Ser['__class__'] = Circuit.__name__
        
        return Ser

    raise TypeError(str(obj) + ' is not JSON serializable')

#Circuit object deserialization function
def deserialize_circuit(obj):
    if '__class__' in obj:
        if obj['__class__'] == 'Circuit':
            Circ = Circuit(obj['name'], obj['miRNA'])
            del obj['__class__']
            Circ.__dict__ = obj
                        
            return Circ

#MAIN
def main(NAME, MIRNA):
    global k, timesuffix, perc_0
        
    #Moment in time:
    timesuffix = '_'.join(
        str(date.now()
        ).split())
    
    ENSEMBLEpre = Circuit(NAME, MIRNA)  
        
    if __name__ == '__main__':
    
        SIMULATION1 = Test_tube(ENSEMBLEpre)
        FUELSIMULATION1 = Test_tube(ENSEMBLEpre, 'FUEL')
    
        OUTFILE = open('Output_{}_{}.txt'
            .format(ENSEMBLEpre.name, timesuffix),
                        'w')
        
        OUTFILE.write('This is the output of your job done on {}\n'
            .format(timesuffix))
        
        OutfileWriter(OUTFILE,
                      ENSEMBLEpre,
                      SIMULATION1,
                      FUELSIMULATION1)

        toolbar = 40
        perc = 0
        print('Progress:')
    
        sys.stdout.write('{}% [{}]\033[92m'
            .format(perc, ' '*toolbar))
    
        sys.stdout.flush()
        sys.stdout.write('\b'*(toolbar+6))

        nums = range(toolbar+1)[1:]
    
    
    #1e5 cycles of mutations and selection following the global score  
    ENSEMBLEpost = Circuit(ENSEMBLEpre.name, ENSEMBLEpre.miRNA)
    ENSEMBLEpost.duplicate(ENSEMBLEpre)
    
    for k in range(int(1+1e5))[1:]: #1e5
        
        ENSEMBLEpost.mutate()
        
        if ENSEMBLEpost.score >= ENSEMBLEpre.score:
            ENSEMBLEpre.duplicate(ENSEMBLEpost)
                    
        else:
            Metropolis(ENSEMBLEpre, ENSEMBLEpost)
            
        if __name__ == '__main__':
            
            perc = int(k*100/(1e5))

            if int(k*toolbar/1e5) in nums:
                sys.stdout.write('{}% [{}{}]'
                    .format(perc,
                            nums[0]*u'\u2588',
                            ' '*(toolbar-nums[0])))
                
                nums = nums[1:]
            
            if k == int(1e5):
                sys.stdout.write('\n\033[0m')

            sys.stdout.flush()
            sys.stdout.write('\b'*(toolbar+6))
    
    SHADOW = shadowcirc(ENSEMBLEpre)
    
    ShadowInfo = ''
    
    for el in SHADOW.__GUIDE__:
          
        setattr(ENSEMBLEpre, el, getattr(SHADOW, el))
        ShadowInfo += '>{}\n{}\n'.format(el, getattr(SHADOW, el))

    try:
        with open('{}_circuit.json'.format(NAME)) as f:
            Circ = json.load(f,
                object_hook=deserialize_circuit)
        
        if ENSEMBLEpre.stdscore > Circ.stdscore:
            with open('{}_circuit.json'.format(NAME), 'w') as f:
                json.dump(ENSEMBLEpre, f,
                    default=serialize_circuit)

    except:
        with open('{}_circuit.json'.format(NAME), 'w') as f:
            json.dump(ENSEMBLEpre, f,
                default=serialize_circuit)
    
    if __name__ == '__main__':    
        SIMULATION1 = Test_tube(ENSEMBLEpre)
        eqsbarplot(SIMULATION1)
    
        FUELSIMULATION1 = Test_tube(ENSEMBLEpre, 'FUEL')
        
        #OUTFILE = open('Output_{}_{}.txt'
         #   .format(ENSEMBLEpre.name, timesuffix),
         #               'w')
        
        #OUTFILE.write('This is the output of your job done on {}\n'
         #   .format(timesuffix))
                
        OutfileWriter(OUTFILE,
                      ENSEMBLEpre,
                      SIMULATION1,
                      FUELSIMULATION1)
            
        OUTFILE.write('\nProposed shadow cancellation circuit\n')
        OUTFILE.write(ShadowInfo)    
        return None

    else:
        return ENSEMBLEpre
        
if __name__ == '__main__':
    #Define initial input
    try:
        USERINPUT = sys.argv[1].upper()
    except:
        USERINPUT = input('Enter your input: '
            ).upper()

    #Checks if input is a raw sequence or a fasta file
    if USERINPUT.lower().split('.')[-1] == 'fasta':
        insequences = fileinput()

        for el in insequences:
            NAME = el
            MIRNA = insequences[el]
            main(NAME, MIRNA)

    else:
        cmdinput()
        main(NAME, MIRNA)

    #NuPACK files cleanup
    TmpCleaner()
        
    elapsed_time = str(
        (time() - start_time)/60)

    print('Job finished on {}. Elapsed time was: {} minutes.'
        .format(str(date.now()).split('.')[0], elapsed_time[:-13]))

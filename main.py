import sys, math
CNTRL = "~~~~~~~"
debug = False


def main(args):
    with open(args[1],'r') as code:
        register, program = sectioner([line for line in code]);
        names,init_state = prepare(register)
        tokens = tokenize(program)
        result = Qal(init_state,tokens,names)
        if debug:
            print('reg: ',register,'\nprog: ',program, '\nnamedict: ',names, '\nprepped: ' , init_state, '\ntokenized: ' ,tokens, '\nresult: ', result)
        

def sectioner(inputfile):
    progFlag, regFlag = False, False
    progFlagUsed, regFlagUsed = False, False
    regout, progout = [],[]
    labels = ['.prog','.reg']

    for line in inputfile:
        
        if regFlag:
           regout.append(line.strip())
        if progFlag:
            progout.append(line.strip())
        
        
        if line.strip() == ".reg":
            if regFlagUsed:
                raise Exception("Register flag used multiple times")
            progFlag, regFlag = False, True
            regFlagUsed = True
            

        if line.strip() == ".prog":
            if progFlagUsed:
                raise Exception("Program flag used multiple times")
            progFlag, regFlag = True, False
            progFlagUsed = True

    return [[i for i in regout if i and i not in labels],[i for i in progout if i and i not in labels]]

def prepare(register):
    #TODO add some better error checking
    nameMapping = {}
    state = []
    for n in range(len(register)):
        line = register[n]
        codes = [i for i in line.split() if i]
        if codes[0][0].lower() not in "qwertyuiopasdfghjklzxcvbnm":
            raise Exception("invalid variable name:\n\n" + line)
        nameMapping[codes[0]] = n

        toKron = []
        if codes[1][0] == "+" and codes[1][1] == "x":
            toKron = [1/2**(1/2),1/2**(1/2)]
        elif codes[1][0] == "-" and codes [1][1] == "x":
            toKron = [-1/2**(1/2),1/2**(1/2)]
        
        elif codes[1][0] == "+" and codes[1][1] == "y":
            toKron = [1/2**(1/2), 0+1j / 2**(1/2)]
        elif codes[1][0] == "-" and codes [1][1] == "y":
            toKron = [1/2**(1/2), 0-1j / 2**(1/2)]
        
        elif codes[1][0] == "+" and codes[1][1] == "z":
            toKron = [1,0]
        elif codes[1][0] == "-" and codes [1][1] == "z":
            toKron = [0,1]

        if len(state) == 0:
            state = toKron
        else:
            state = kron_vector(state,toKron) 
    return [nameMapping,state]

def kron_vector(a,b):
    out = []
    for alpha in a:
        out+=[alpha*beta for beta in b]
    return out

def kron_matrix(m1,m2):
    r1,c1 = len(m1), len(m1[0])
    r2,c2 = len(m2), len(m2[0])
    
    def prod(a1,a2,i1,j1,i2,j2):
        if a1 is CNTRL:
            return CNTRL if i2==j2 else 0
        if a2 is CNTRL:
            return CNTRL if i1==j1 else 0
        else:
            return a1*a2
    return [[prod(m1[i1][j1],m2[i2][j2],i1,j1,i2,j2) for i1 in range(r1) for i2 in range(r2)] for j1 in range(c1) for j2 in range(c2)] #im sorry


def tokenize(prog):
    #TODO add some better error checking
    out = []
    for line in prog:
        codes = [i for i in line.split(' ',1) if i]
        instr = codes[0]
        try:
            arguments = codes[1].replace(' ', '').split(',')
        except :
            arguments = []
        out += [[instr,arguments]]
    return out






def Qal(state, tokens, names):
    lowercase_instr_set = {'h':1,'cx':2,'s':1,'t':1, 'st':0}
    onearg_instr_to_matrix = {'h': [[1/math.sqrt(2),1/math.sqrt(2)],[1/math.sqrt(2),-1/math.sqrt(2)]],
                              's': [[1,0],[0,0+1j]],
                              't': [[1,0],[0,math.e**(0+1j*math.pi/4)]],
                              'x': [[0,1],[1,0]],
                              'y': [[0,0-1j],[0+1j,0]],
                              'z': [[1,0],[0,-1]]}

    names = dict(names)
    for line in tokens:
        instr = line[0]
        argus = line[1]
        
        if instr.lower() not in lowercase_instr_set:
            raise Exception("invalid instruction:\n" + line)
        if len(argus) != lowercase_instr_set[instr.lower()]:
            raise Exception("incorrect amount of arguments:\n"+ line)
        
        if instr.lower() in ['h','s','t','x','y','z']:
            ind = names[argus[0]]
            gate = onearg_instr_to_matrix[instr.lower()]

            toMult = gate if ind == 0 else [[1,0],[0,1]]   #set the first part of the unitary up
            for i in range(1,int(math.log(len(state),2))): #for i in range(1,numqubits)
                if i == ind:
                    toMult = kron_matrix(toMult,gate)      #compute kronecker product of our gate for the target qubit's 
                else:
                    toMult = kron_matrix(toMult, [[1,0],[0,1]])
            
            state = [sum(q) for q in [[state[j] * toMult[i][j] for j in range(len(toMult[0]))] for i in range(len(toMult))]] #sorry again. computes U|psi>
        
        elif instr.lower() == 'cx':
            C = [[CNTRL,0],[0,1]]
            X = [[0,1],[1,0]]
            I = [[1,0],[0,1]]
            ctrlqubit = names[argus[0]]
            targetqubit  = names[argus[1]]
            toMult = X if targetqubit == 0 else (C if ctrlqubit == 0 else I)
            for i in range(1,int(math.log(len(state),2))):
                if i == ctrlqubit:
                    toMult = kron_matrix(toMult,C)
                elif i == targetqubit:
                    toMult = kron_matrix(toMult,X)
                else:
                    toMult = kron_matrix(toMult, [[1,0],[0,1]])
            toMult = [[(i if i is not CNTRL else 1) for i in j] for j in toMult]
           
            state = [sum(q) for q in [[state[j] * toMult[i][j] for j in range(len(toMult[0]))] for i in range(len(toMult))]] #sorry again. computes U|psi>
        
        #output instructions
        elif instr.lower() == 'st':
            print('[INFO] the "st" instruction has been called:\n',state)


        else:
            raise Exception("invalid instruction?? how did you get here")
    return state
            



main(sys.argv)

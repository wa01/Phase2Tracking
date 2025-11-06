#
# Definitions for histograms from a configuration file
#
from fnmatch import fnmatch

class HistogramDefinition:
    ''' A single histogram definition. 
    '''
    reqGenFields = [ 'canvasName', 'histogramName', 'histogramTitle', 'variable', 'baseCuts' ]
    reqHistFields = [ ]
    requiredFields = reqGenFields + reqHistFields
    optGenFields =  [ 'effCuts', 'logY' ]
    optHistFields = [ 'xNbins', 'xMin', 'xMax', 'xTitle', 'yTitle', 'yNbins', 'yMin', 'yMax', \
                          'zMin', 'zMax', 'display', 'profile' ]
    optionalFields = optGenFields + optHistFields
    allFields = requiredFields + optionalFields
    allHistFields = reqHistFields + optHistFields

    def __init__(self,name,inputDict):
        ''' Define histogram and drawing parameters from dictionary.
            Arguments:
              name ....... name of the histogram definition
              inputDict .. dictionary with values
        '''
        #
        assert name.isalnum()
        self.name = name
        #
        # read parameters
        #
        # default is None for all parameters
        #
        self.parameters = { x:None for x in HistogramDefinition.allFields }
        #
        # loop over main dictionary
        #
        for k,v in inputDict.items():
            #
            # skip standard python entries
            #
            if k.startswith('__'):
                continue
            #
            # standard entry in main dictionary - store value
            #
            if k in HistogramDefinition.allFields:
                #
                # general variable
                #
                self.parameters[k] = v
            #
            # nested dictionary with parameters for a specific module type
            #
            elif k.startswith("mType"):
                #
                # mType-specific histogram parameters
                #
                assert type(v)==dict and ( not k in self.parameters ) and len(k)>5 and k[5:].isdigit()
                self.parameters[k] = { x:None for x in HistogramDefinition.allHistFields }
                for kh,vh in v.items():
                    if kh in HistogramDefinition.allHistFields:
                        self.parameters[k][kh] = vh
                    else:
                        print("Warning: key",kh, \
                                  "is not a standard histogram field name - ignoring the entry in", \
                                  self.name)
            #
            # unknown key
            #
            else:
                print("Warning: key",k,"is not a standard field name - ignoring the entry in",self.name)
        #
        # set some parameters from other inputs if unspecified
        #
        if self.parameters['canvasName']==None:
            self.parameters['canvasName'] = "c" + self.name[0].upper() + self.name[1:]
        if self.parameters['histogramName']==None:
            self.parameters['histogramName'] = "h" + self.name[0].upper() + self.name[1:]
        #
        # make sure all required fields are present
        #
        for f in HistogramDefinition.requiredFields:
            assert ( f in self.parameters ) and self.parameters[f]!=None

    #def __getitem__(self,field):
    #    if field in self.parameters:
    #        return self.parameters[field]
    #    return None

    def getParameter(self,name,mType=None):
        ''' Retrieve a single parameter (optionally a specific one for a given module type)
        '''
        result = None
        #
        # give priority to parameter specific to a module type
        #
        mTName = "mType"+str(mType) if mType!=None else None
        if ( mTName in self.parameters ) and ( name in self.parameters[mTName] ):
          result = self.parameters[mTName][name]
          if result!=None:
            return self.parameters[mTName][name]
        #
        # not found: use general parameter
        #
        if name in self.parameters:
            return self.parameters[name]
        return None

    def vetoMType(self,mType):
        ''' Check for an mType entry with display = False
        '''
        #
        # try to get 'display' parameter
        #
        name = 'display'
        mTName = "mType"+str(mType) if mType!=None else None
        if ( mTName in self.parameters ) and ( name in self.parameters[mTName] ):
            if self.parameters[mTName][name]!=None:
                return not self.parameters[mTName][name]
        return False

class HistogramDefinitions:

    def __init__(self):
        self.allDefinitions = { }
        self.allHistoNames = set()
        self.allCanvases = set()

    def add(self,hdef):
        assert not hdef.name in self.allDefinitions
        assert not hdef.getParameter('histogramName') in self.allHistoNames
        assert not hdef.getParameter('canvasName') in self.allCanvases
        self.allDefinitions[hdef.name] = hdef
        self.allHistoNames.add(hdef.getParameter('histogramName'))
        self.allCanvases.add(hdef.getParameter('canvasName'))
        

    def __getitem__(self,name):
        if name in self.allDefinitions:
            return self.allDefinitions[name]
        return None
        
def loadHistogramDefinitions(configName,selectedNames=[],vetoedNames=[]):
    ''' Load histogram definitions from a configuration file. The configuration file
        defines dictionaries with the definitions as variable. The name of the variable
        serves as name of the HistogramDefinition.
        Arguments:
           configName ..... name of the module to import from / python file name
           selectedNames .. explicit list of histogram names to be imported (filename wildcard syntax)
           vetoedNames .... explicit list of histogram names to be skipped (filename wildcard syntax)
    '''
    # load histogram definitions
    #
    result = HistogramDefinitions()
    if configName==None:
        return result

    moduleName = configName[:-3] if configName.endswith(".py") else configName
    module = __import__(moduleName)
    for n in dir(module):
        if n.startswith('__'):
            continue
        hDict = getattr(module,n)
        #print(n,type(hDict))
        #sys.exit()
        assert type(hDict)==dict
        #
        # check if in list of histograms to be displayed
        #
        selFlg = False
        for p in selectedNames:
            if fnmatch(n,p):
                selFlg = True
                break
        if not selFlg:
            continue
        #
        # check if in list of histograms to be vetoed
        #
        selFlg = True
        for p in vetoedNames:
            if fnmatch(n,p):
                selFlg = False
                break
        if not selFlg:
            continue
        #
        # add histogram
        #
        hDef = HistogramDefinition(n,hDict)
        result.add(hDef)
        print("Added",hDef.getParameter('canvasName'))

    return result


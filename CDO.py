import bezier
import numpy as np
import subprocess
import os
import shutil
import collections
import pandas as pd
import csv
from pyDOE import lhs
import matplotlib.pyplot as plt

# Class for computational design optimisation study of a centrifugal pump
class CDO:
#   Initialise some important class objects
    def __init__(self,Fixed=None,Nominal=None,MeshParam=None,CFDParam=None):
#       The directory in which the study will be performed in
#       NOTE: This directory must contain template files: TempCheckMesh.rpl,
#       TempImpellerCFD.ccl, TempImpellerGeom.bgi, TempImpellerMesh.tse, and
#       TempImpellerSetup.pre.
        self.Dir = 'C:\Users\james\Documents\BladeGen'
#       The Ansys directory containing all of the executable and bat files
#       necessary to call the applications to create, mesh, and evaluate the
#       geometries i.e. BladeBatch.exe, cfxtg.exe, cfx5pre.exe, icemcfd.bat,
#       and cfx5solve.exe.
        self.AnsysDir = r'C:\Program Files\ANSYS Inc\v182'
#       The dictionary of geometry parameters that will be fixed throughout the
#       study e.g. For pump impeller, inlet and impeller radii will often be 
#       fixed.
        if Fixed is not None:
            self.Fixed = Fixed
        else:
            self.Fixed = {'R1': 50.988300, 'R2': 97.622400, 'Rcut': 12.065000, \
         'Hin': 100.0, 'Rout': 120.0, 'HBz3': 0.0, 'HBz4': 0.0, 'SBz4': 22.569800,\
         'SBr2': 50.988300,'bladeNo': 6} 
#       The dictionary of geometry parameters for a nominal geometry. The value
#       of these parameters will change from geometry to geometry, however this
#       dictionary will be left untouched.
        if Nominal is not None:
            self.Nominal = Nominal
        else:
            self.Nominal = collections.OrderedDict()
            self.Nominal['Hlenfrac'] = 0.461192
            self.Nominal['Slenfrac'] = 0.359983
            self.Nominal['SBz1'] = 40.861900
            self.Nominal['SBz2'] = 32.602400
            self.Nominal['SBz3'] = 26.490900
            self.Nominal['SBr3'] = 56.544300
            self.Nominal['HBz1'] = 22.472400
            self.Nominal['HBz2'] = 8.225470
            self.Nominal['HBr2'] = 20.290500
            self.Nominal['HBr3'] = 34.537400
            self.Nominal['LEz2'] = 0.0
            self.Nominal['LEr2'] = 0.0
            self.Nominal['LEz3'] = 0.0
            self.Nominal['LEr3'] = 0.0
            self.Nominal['Beta11'] = 34.831
            self.Nominal['Beta12'] = 50.357
            self.Nominal['Beta13'] = 45.441
            self.Nominal['Beta14'] = 22.500
            self.Nominal['Beta51'] = 21.018
            self.Nominal['Beta52'] = 22.500
            self.Nominal['Beta53'] = 22.994
            self.Nominal['Beta54'] = 22.500
            self.Nominal['Thk11'] = 5.857
            self.Nominal['Thk12'] = 5.857
            self.Nominal['Thk13'] = 5.857
            self.Nominal['Thk14'] = 5.857
            self.Nominal['Thk51'] = 5.857
            self.Nominal['Thk52'] = 5.857
            self.Nominal['Thk53'] = 5.857
            self.Nominal['Thk54'] = 5.857
#       The dictionary of meshing parameters for the TSE file used in TurboGrid
        if MeshParam is not None:
            self.MeshParam = MeshParam
        else:
            self.MeshParam = {'bladeNo': self.Fixed['bladeNo'], 'PeriodicBias': 0.6, 
             'InHubLoc': 0.1, 'InShroudLoc': 0.2, 'OutletLoc': 0.5, \
             'FiveVertFac': 1.0, 'GlobSizeFac': 1.0, 'BLFacBase': 3.0, \
             'BLFacRate': 0.0, 'MaxExpRate': 1.3, 'InGrow': 1.05, \
             'OutGrow': 1.03, 'SpanFac': 1.08, 'TESplitFac': 1.0}
#       The dictionary of CFD parameters for the CCL file used for CFX Pre
        if CFDParam is not None:
            self.CFDParam = CFDParam
        else:
            self.CFDParam = {'AngularVel': 261.8, 'FlowRateKGPerS': 30.0, \
            'PressureBC': 600000,  'MaxIter': 800, 'noPartitions': 8, \
            'ConvergenceRes': 0.00001}

#   Creates a list of 50 spline points from a list of Bezier control points 
#   defining a 2D line using the library "bezier"
#   Inputs:
#   Z:          a list of Bezier CP values for the first dimension
#   R:          a list of Bezier CP values for the second dimension
#   lenfrac:    a float with values 0-1 that defines how far along curve
#               the leading edge will begin
#   Outputs:
#   Points:     a list of 50 spline points defining the curve
#   LE:         a 2D array specifying the location of the LE
    @staticmethod
    def Bezier2SplinePoints(Z,R,lenfrac):
#       Number of desired spline points minus one (does not include LE point)
    	n = 49
#       Define the nodes as a Fortran array in order to use "bezier" library
    	nodes = np.asfortranarray([Z,R,])
#       Define the Bezier curve from the nodes (degree has been set to 3 as
#       we have specified 4 CPs)
    	curve = bezier.Curve(nodes, degree=3)
#       Define percentages along the curve to take spline points from, includes
#       the length fraction
    	values = np.sort(np.append(np.linspace(0,1,n),lenfrac))
#       Initalise the list of spline points and iterate through values to find
#       the points at the specified percentages of the curve
    	Points = []
    	for ind,val in enumerate(values):
    		point = (curve.evaluate(val))
    		Points.append(point)
#           Define the point in which the percentage is equal to the given 
#           length fraction i.e. the leading edge point
    		if val == lenfrac:
    			LE = point
    	return Points,LE
    
#   Modifies a file by updating all the strings in a file matching the strings
#   in a given dictionary with the value associated with the key
#   Inputs:
#   FileIn:     Input file
#   FileOut:    Output file
#   Dict:       Dictionary with key values being the string to be swapped with 
#               values associated to that key
    @staticmethod
    def UpdateFileFromDict(FileIn,FileOut,Dict):
#       Open input file and read data
        with open(FileIn,'r') as f:
            filedata = f.read()
#       Remove the input file once data has been read in
        os.remove(FileIn)
#       Initialise the new file data and iterate through keys and values in 
#       dictionary, replacing all keys found in file data with the value
#       associated with the key
        newdata = filedata
        for key,value in Dict.items():
            newdata = newdata.replace(key,str(value))
#       Open output file and write the new data
        with open(FileOut,'w') as f:
            f.write(newdata)
        return None
    
#   Takes a list of dictionaries with all the same key values and writes a CSV
#   with the header equal to the key values and each containing the values
#   for each dictionary instance
#   Inputs:
#   FileOut:    Output file
#   Dicts:      List of dictionaries
    @staticmethod
    def DictList2CSV(FileOut,Dicts):
#       Find the key values of the first dictionary 
#       NOTE: all dictionaries must have the same keys!
        keys = Dicts[0].keys()
#       Open the output file and write the entries of each dictionary as a row
#       in the CSV
        with open(FileOut, 'wb') as outputfile:
            dictwriter = csv.DictWriter(outputfile, keys)
            dictwriter.writeheader()
            dictwriter.writerows(Dicts)
        return None

#   Appends an output file with results from an ordered dictionary
#   NOTE: The input must be an ordered dictionary, as the values are inserted
#   in the order they are placed originally
#   Inputs:
#   FileOut:    Output file, which must be a CSV
#   Dict:       Ordered dictionary with values corresponding to the same keys
#               as any previous dictionaries that have been used
    @staticmethod
    def OrderedDict2CSV(FileOut,Dict):
#       Check if the file exists
        if os.path.exists(FileOut):
            tag = False
        else:
            tag = True
#       Open file and append results. If file does not exist, both keys and 
#       values are input, whereas if it does exist, only values are input.
        with open(FileOut,'ab') as f:
            w = csv.writer(f)
            if tag:
                w.writerow(Dict.keys())
            w.writerow(Dict.values())
        return None
    
#   Extracts the results with header defined in a list from a CSV file using the
#   "pandas" library
#   Inputs:
#   FileName:   Name of CSV file
#   Tags:       List of header tags to extract data from
#   Outputs:
#   Results:    Dictionary of results with header values contained in Tags
    @staticmethod
    def ExtractResultsFromCSV(FileName,Tags):
#       Read in CSV file
        df = pd.read_csv(FileName)
#       Initialise dictionary and update the values with the columns from the 
#       CSV file
        Results = {}
        for num,key in enumerate(Tags):
            Results[key] = df[key].values
        return Results
    
#   'Cleans' a data set read from relational data table stored in CSV format
#   by removing any row with NaN values and normalising the columns
#   Inputs:
#   FileName:   Name of CSV file
#   Outputs:
#   df:         Dataframe object containing all the cleaned data
    @staticmethod
    def CleanData(FileName,scale=[0,1]):
        df = pd.read_csv(FileName)
        df = df.dropna(thresh=len(df.columns))
#        df[df.columns[:30]] = scale[0]+(scale[1]-scale[0])*(df[df.columns[:30]]-df[df.columns[:30]].min())/(df[df.columns[:30]].max()-df[df.columns[:30]].min())
        df=scale[0]+(scale[1]-scale[0])*(df-df.min())/(df.max()-df.min())
        return df
    
    @staticmethod
    def CSVDataFrame2Dict(FileName):
        Data = pd.read_csv(FileName)
        return Data.to_dict('records')[0]
    
    def CreateWorkDir(self,FileName):
        self.WorkDir = self.Dir+'/'+FileName
        if os.path.exists(self.WorkDir):
            shutil.rmtree(self.WorkDir)
        os.mkdir(self.WorkDir)
        return None
    
#   Creates interactive parallel coordinates given a relational data table in  
#   CSV form using the library "plotly"
#   NOTE: a credentials file must located in the "home" directory with username
#   and API key information to use this library
#   E.g. username = james.gross189 api_key = ksQEl0SV0y3qmIZXfpJ1
#   Inputs:
#   FileName:   the location of the relational data table in CSV form
    def ParallelCoord(self,FileName):
#       Import plotly library (set credentials if necessary)
#        plotly.tools.set_credentials_file(username=ID, api_key=key)
        import plotly.plotly as py
        import plotly.graph_objs as go
#       Read in the data from the CSV file, drop rows with NaN values in any 
#       entry, and normalise the entries
        df = self.CleanData(FileName)
#        df = pd.read_csv(FileName)
#        df = df.dropna(thresh=len(df.columns))
#       Initialise the data and loop through the columns in the dataframe 
#       creating a dictionary containing necessary info for plotly parallel 
#       coordinates
        Data = []
        for column in df:
            Data.append(dict(range=[0,1],label=column, \
                             values = df[column])) 
        data = [go.Parcoords(line = dict(color = '#25fc25'),dimensions = list(Data))]
#       Plot the coordinates and open up a webpage to view the interactive
#       coordinate system
        fig = go.Figure(data = data)
        py.iplot(fig, filename = 'ParallelCoords', auto_open=True)
        return None
#    
#   Creates the dictionary used when updating the BGI file which is used by
#   the "BladeBatch.exe" i.e. write input file for BladeModeler
#   Inputs:
#   Param:      Dictionary containing all parameter values defining the 
#               geometries 
#               NOTE: The keys should be the same as that as the self.Nominal
#               dictionary, but with possibly different values
    def WriteBGI(self,Param):
#       Define input and output files and initialise the dictionary used to 
#       modify the BGI file
        FileIn = self.WorkDir+'\TempImpellerGeom.bgi'
        FileOut = self.WorkDir+'\ImpellerGeom.bgi'
        Dict = {}
#       Add all entries from the self.Fixed dictionary
        for key,value in self.Fixed.items():
            Dict[key] = value
#       Add all entries from the Param dictionary
        for key,value in Param.items():
            Dict[key] = value
#       Create a list of R and Z values of the Bezier CPs defining both the 
#       hub and shroud curves
        HZB = [Param['HBz1'], Param['HBz2'], self.Fixed['HBz3'], self.Fixed['HBz4']]
        HRB = [self.Fixed['Rcut'], Param['HBr2'], Param['HBr3'], self.Fixed['R2']]
        SZB = [Param['SBz1'], Param['SBz2'], Param['SBz3'], self.Fixed['SBz4']]
        SRB = [self.Fixed['R1'], self.Fixed['SBr2'], Param['SBr3'], self.Fixed['R2']]
#       Find both the hub and shroud spline points, as well as the leading
#       edge points
        HPoints,LE1 = self.Bezier2SplinePoints(HZB,HRB,Param['Hlenfrac'])
        SPoints,LE4 = self.Bezier2SplinePoints(SZB,SRB,Param['Slenfrac'])
        delLE = 0.1*np.linalg.norm(LE4-LE1)
#       Add all of the spline points for both curves to the dictionary
        for num,point in enumerate(HPoints):
            Dict['H'+str(num+1)+'z'] = point[0,0]
            Dict['H'+str(num+1)+'r'] = point[1,0]
        for num,point in enumerate(SPoints):
            Dict['S'+str(num+1)+'z'] = point[0,0]
            Dict['S'+str(num+1)+'r'] = point[1,0]
#       Calculate points for LE2 and LE3
        LE2 = LE1 + (LE4-LE1)/3.
        LE3 = LE1 + 2.*(LE4-LE1)/3.
#       Add the leading edge points to the dictionary
        Dict['LEz1'] = LE1[0,0]
        Dict['LEr1'] = LE1[1,0]
        Dict['LEz2'] = LE2[0,0] + Param['LEz2']*delLE
        Dict['LEr2'] = LE2[1,0] + Param['LEr2']*delLE
        Dict['LEz3'] = LE3[0,0] + Param['LEz3']*delLE
        Dict['LEr3'] = LE3[1,0] + Param['LEr3']*delLE
        Dict['LEz4'] = LE4[0,0]
        Dict['LEr4'] = LE4[1,0] 
#       Modify the template BGI file with dictionary values
        self.UpdateFileFromDict(FileIn,FileOut,Dict)
        return None

#   Creates the dictionary used when updating the TSE file which is used by
#   the "cfxtg.exe" i.e. write input file for TurboGrid
    def WriteTSE(self):
#       Define input and output files and initialise the dictionary used to 
#       modify the TSE file
        FileIn = self.WorkDir+'\TempImpellerMesh.tse'
        FileOut = self.WorkDir+'\ImpellerMesh.tse' 
#       Creates dictionary and initialises with values from the self.MeshParam 
#       dictionary
        Dict = self.MeshParam.copy()
#       Add the file names for all of the curves and the mesh to the dictionary
        Dict['HubFile'] = self.WorkDir+'\hub.curve'
        Dict['ShroudFile'] = self.WorkDir+'\shroud.curve'
        Dict['BladeFile'] = self.WorkDir+'\profile.curve'
        Dict['MeshFile'] = self.WorkDir+'\Mesh.gtm'
#       Modify the template TSE file with dictionary values        
        self.UpdateFileFromDict(FileIn,FileOut,Dict)
        return None

#   Creates the dictionary used when updating the PRE file which is used by
#   the "cfx5pre.exe" i.e. write first input file for CFX Pre
    def WritePRE(self):
#       Define input and output files and initialise the dictionary used to 
#       modify the PRE file
        FileIn = self.WorkDir+'\TempImpellerSetup.pre'
        FileOut = self.WorkDir+'\ImpellerSetup.pre'
        Dict = {}
#       Add the file names for mesh, CCL, and DEF files to the dictionary       
        Dict['MeshFile'] = self.WorkDir+'\Mesh.gtm'
        Dict['CCLFile'] = self.WorkDir+'\ImpellerCFD.ccl'
        Dict['DEFFile'] = self.WorkDir+'\ImpellerCFD.def'
#       Modify the template PRE file with dictionary values         
        self.UpdateFileFromDict(FileIn,FileOut,Dict)
        return None

#   Creates the dictionary used when updating the CCL file which is used by
#   the "cfx5pre.exe" i.e. write second input file for CFX Pre    
    def WriteCCL(self):
#       Define input and output files and initialise the dictionary used to 
#       modify the CCL file
        FileIn = self.WorkDir+'\TempImpellerCFD.ccl'
        FileOut = self.WorkDir+'\ImpellerCFD.ccl'
#       Creates dictionary and initialises with values from the self.CFDParam 
#       dictionary
        Dict = self.CFDParam.copy()
#       Add the flow rate used in the head calculation in kg/s (for water at 25C)
#       NOTE: This is divided by the total blade number 
        Dict['FlowRate4Calc'] = Dict['FlowRateKGPerS']/(self.Fixed['bladeNo']*997.0)
#       Modify the template CCL file with dictionary values   
        self.UpdateFileFromDict(FileIn,FileOut,Dict)
        return None

#   Creates the dictionary used when updating the RPL file which is used by
#   the "icemcfd.bat" i.e. write input file for ICEMCFD    
    def WriteRPL(self):
#       Define input and output files and initialise the dictionary used to 
#       modify the RPL file
        FileIn = self.WorkDir+'\TempCheckMesh.rpl'
        FileOut = self.WorkDir+'\CheckMesh.rpl'
        Dict = {}
#       Define the DEF file and report file name and location
        DEFFileTemp =  self.WorkDir+'\ImpellerCFD.def'
        ReportFileTemp = self.WorkDir+'\MeshReport.txt'
#       Swap all backslashes with forward slashes (necessary for ICEMCFD scripting)
        Dict['DEFFile'] = DEFFileTemp.replace('\\', '/')
        Dict['ReportFile'] = ReportFileTemp.replace('\\', '/')
#       Modify the template RPL file with dictionary values 
        self.UpdateFileFromDict(FileIn,FileOut,Dict)
        return None
    
#   Checks the mesh report against 3 mesh quality metrics i.e. cell quality
#   (defined by ICEM CFD), minimum determinant value, and minimum angle values
#   Outputs:
#   True or False:  True if mesh quality is adequate, False if it is not
#   TODO: Extract number of elements too!! Check out BeautifulSoup?
    def CheckMesh(self):
#       Open mesh report file and reads in data
        FileName = self.WorkDir+'\MeshReport.txt'
        with open(FileName,'r') as f:
            lines = f.readlines()
        for num,line in enumerate(lines):
#           Find the table for cell quality values and the minimum 
#           cell quality and checks it against arbitrary threshold value
            if 'Quality values' in line:
                minQual = float(lines[num+1].split(' ')[1])
                if minQual < 0.18:
                    return False
#           Find the table for determinant values and the minimum 
#           determinant and checks it against arbitrary threshold value
            elif 'Determinant values' in line:
                minDeter = float(lines[num+1].split(' ')[1])
                if minDeter < 0.2:
                    return False
#           Find the table for minimum angle values and the absolute  
#           minimum angle and checks it against arbitrary threshold value
            elif 'Min angle values' in line:
                minAngle = float(lines[num+1].split(' ')[1])
                if minAngle < 10.0:
                    return False
            else:
                continue
#       Return True if mesh is adequate
        return True
        
    def ExtractMeshSize(self):
        FileName = self.WorkDir+'\MeshReport.txt'
        if os.path.exists(FileName):
            with open(FileName,'r') as f:
                lines = f.readlines()
            return int(lines[13].split(' ')[1])
        else:
            return -99999
            
#   Perturbs the input values found in the self.Nominal dictionary using a 
#   given maximum perturbation and perturbation value between 0 and 1
#   Inputs:
#   Perturb:    Dictionary of perturbation values for each input (must be between
#               0 and 1)
#   MaxPerturb: (optional) Dictionary containing all the maximum allowed 
#               perturbation values for each value in self.Nominal dictionary
#   Outputs:
#   Param:      Dictionary containing all perturbed parameter values 
    def PerturbParam(self,Perturb,MaxPerturb):
#       Initialise Param dictionary and update values by perturbing the nominal
#       value such that the new value lies within a circle centred around the 
#       nominal with radius equal to MaxPerturb
        Param = collections.OrderedDict()
        for key,val in self.Nominal.items():
            Param[key] = Perturb[key]*MaxPerturb[key]+val
        return Param

#   Create and evaluate the geometry defined by the dictionary Param using 
#   Ansys packages (i.e. calls BladeModeler, TurboGrid, ICEM CFD, and CFX) by 
#   creating a working directory, copying templates and modifying the entries
#   in the new directory, and then sequentially calling the necessary commands
#   in Ansys
#   Inputs:
#   Param:  Dictionary containing all the parameter values allowed to vary in
#           subsequent runs
#   num:    An integer which is used to create the new working directory 
    def Evaluate(self,Param):
#       Define tags for quantities of interest
        QOI = ['H','nu']
        Results = collections.OrderedDict()
#       Copy all of the template files from starting directory to the working
#       directory
        for FileName in os.listdir(self.Dir):
            if FileName.startswith('Temp'):
                shutil.copy(self.Dir+'\\'+FileName,self.WorkDir)
#       Write the BGI file to the working directory and call the BladeModeler
#       command
        self.WriteBGI(Param) 
        subprocess.call([self.AnsysDir+r"\aisol\BladeModeler\BladeGen\BladeBatch.exe", \
                         self.WorkDir+'\ImpellerGeom.bgi', '-TG', self.WorkDir])
#       Write the TSE file to the working directory and call the TurboGrid
#       command        
        self.WriteTSE()
        subprocess.call([self.AnsysDir+r"\TurboGrid\bin\cfxtg.exe", '-batch', \
                         self.WorkDir+'\ImpellerMesh.tse'])
#       Write the PRE and CCL files to be used for CFX Pre
        self.WritePRE()
        self.WriteCCL()
#       Call CFX Pre and handle call errors by returning nan values for results
        try:
            subprocess.call([self.AnsysDir+r"\CFX\bin\cfx5pre.exe", '-batch', \
                             self.WorkDir+'\ImpellerSetup.pre'])
        except subprocess.CalledProcessError:
            for key in QOI:
                Results[key]=float('nan')
            return Results
#       Write RPL file to working directory and call ICEM CFD to check mesh
#       quality
        self.WriteRPL()
        subprocess.call([self.AnsysDir+r"\icemcfd\win64_amd\bin\icemcfd.bat", \
                         '-batch', self.WorkDir+'\CheckMesh.rpl'])
#       If mesh does not pass mesh quality check, return nan values
        if not self.CheckMesh():
            for key in QOI:
                Results[key]=float('nan')
#       If mesh passes mesh quality check, continue evaluation procedure
        else:
#           Call CFX solver in working directory using the number of partitions
#           defined in self.CFDParam
            subprocess.call([self.AnsysDir+r"\CFX\bin\cfx5solve.exe", '-chdir', \
                             self.WorkDir, '-def', self.WorkDir+'\ImpellerCFD.def', \
                             '-par-local', '-partition', \
                             str(self.CFDParam['noPartitions'])])
#           Handle call errors when attempting to extract results by returning
#           nan values 
#           NOTE: call errors typically happen here when CFX5SOLVE could run,
#           however did not return any results (e.g. when cell elements have
#           negative volumes)
            try:
#               Temporarily save the results to workspace so that errors may
#               effectively be handled
                Monitor = subprocess.check_output([self.AnsysDir+r"\CFX\bin\cfx5mondata.exe",\
                                         '-res', self.WorkDir+'\ImpellerCFD_001.res'])
            except subprocess.CalledProcessError:
                for key in QOI:
                    Results[key]=float('nan')
                return Results
#           Save the monitor results to a CSV file to extract the necessary
#           data for inspection
            Filename = self.WorkDir+'\Monitor.csv'
            with open(Filename,'w') as f:
                f.write(Monitor)
#           Define the column tags in CSV file for the quantities of interest 
#           NOTE: This is because cfx5mondata stores monitor values as 
#           "USER POINT, varName"
            Columns = []
            for item in QOI:
                Columns.append('USER POINT,'+item)
#           Extract the results from the CSV file
            Res = self.ExtractResultsFromCSV(Filename,Columns)
#           Update the Param dictionary with the extracted results from the most
#           recent timestep
            for key,val in Res.items():
                ID = key.split(',')[-1]
                if len(val)<self.CFDParam['MaxIter']:
                    Results[ID]=val[-1]
                else:
                    Results[ID]=float('nan')
        return Results
    
#   Runs a DOE study by taking a LHS of the design variables over a specified
#   number of sample and iteratively perturbing parameters around these values
#   to acquire new geometries, which are then evaluated using the "Evaluate"
#   method
#   Inputs:
#   noSamples:  number of samples to run DOE study for
#   NOTE: This value must be an integer greater than 1
    def DOE(self,noSamples,lhd=None):
        if lhd is not None:
            Samples = []
            for i in range(len(lhd)):
                Sample = collections.OrderedDict()
                for num,val in enumerate(self.Nominal):
                    Sample[val] = lhd[i,num]
                Samples.append(Sample)
            MaxPerturb = self.CSVDataFrame2Dict(self.Dir+'/MaxPerturb.csv')
        else:
#           Latin Hypercube samples
            lhd = -1.+2.*lhs(len(self.Nominal), samples=noSamples)
            Samples = []
            for i in range(len(lhd)):
                Sample = collections.OrderedDict()
                for num,val in enumerate(self.Nominal):
                    Sample[val] = lhd[i,num]
                Samples.append(Sample)
#           Define maximum perturbations for each input
            MaxPerturb = self.Nominal.copy()
            for key,val in MaxPerturb.items():
                MaxPerturb[key] = 0.1*val
            MaxPerturb['LEz2'] = 1.0
            MaxPerturb['LEr2'] = 1.0
            MaxPerturb['LEz3'] = 1.0
            MaxPerturb['LEr3'] = 1.0
#           Write max perturbation values to file
            self.DictList2CSV(self.Dir+'\MaxPerturb.csv',[MaxPerturb])
#           Write nominal geometry values to file
            self.DictList2CSV(self.Dir+'\Nominal.csv',[self.Nominal])
#           Write all the LHS points to a file
            self.DictList2CSV(self.Dir+'\Sample.csv',Samples)
#       Loop over the samples, perturbing the nominal values, and evaluating
#       the geometries that are returned
        for i in range(len(lhd)):
            Param = self.PerturbParam(Samples[i],MaxPerturb)
#           Define and create new working directory if it does not exist already.
#           If it exists already, delete it and make a new directory
            self.CreateWorkDir('Work'+str(i))
            Results = self.Evaluate(Param)
            Output = Samples[i].copy()
            for key,value in Results.items():
                Output[key] = value
            shutil.rmtree(self.WorkDir)
#           Print the results to a CSV file
            self.OrderedDict2CSV(self.Dir+'\Results.csv',Output)
        return None
    
#   Performs a grid independence study of the nominal geometry by updating the
#   TurboGrid mesh parameter 'Global Size Factor' and saves the results in a
#   CSV file for later investigation
    def GridIndependenceStudy(self):
#       Define the range of values to loop over for the grid independence study
#        GlobSizeFac = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
#        GlobSizeFac = [1.4,1.5]
#        for num, val in enumerate(GlobSizeFac):
##           Initialise the Results and Param dictionaries
#            Results = {}
#            Param = self.Nominal.copy()
##           Update the Global Size Factor mesh parameters in the self.MeshParam
#            self.MeshParam['GlobSizeFac'] = val
##           Define and create new working directory if it does not exist already.
##           If it exists already, delete it and make a new directory  
#            self.WorkDir = self.Dir+'\Study'+str(num)
#            if os.path.exists(self.WorkDir):
#                shutil.rmtree(self.WorkDir)
#            os.mkdir(self.WorkDir) 
#            Param = self.Evaluate(Param)
##           Update the Results dictionary with the CFD results
#            Results['GlobSizeFac'] = val
#            Results['NoElements'] = self.ExtractMeshSize()
#            Results['H'] = Param['H']
#            Results['nu'] = Param['nu']
##           Write results to file
#            self.OrderedDict2CSV(self.Dir+'\GridIndependence.csv',Results)
##       Collect the data from the CSV file and plot the results from study
        Results = self.ExtractResultsFromCSV(self.Dir+'\GridIndependence.csv',['NoElements','H','nu'])
        plt.scatter(Results['NoElements']/1000,Results['nu'])
        plt.show()
        return None 

    def EvaluateSample(self,Array,num):
        self.CreateWorkDir('New'+str(num))
        Sample = collections.OrderedDict()
        for num,val in enumerate(self.Nominal):
            Sample[val] = Array[num]
        MaxPerturb = self.CSVDataFrame2Dict(self.Dir+'/MaxPerturb.csv')
        Param = self.PerturbParam(Sample,MaxPerturb)
        Results = self.Evaluate(Param)
        Output = Sample.copy()
        for key,value in Results.items():
            Output[key] = value  
        self.OrderedDict2CSV(self.Dir+'\AddResults.csv',Output)
#        print Output
        return None
            
if __name__ == '__main__':
    ins = CDO()
#    ins.EvaluateSample(x)
############################### Grid Study ##################################
#    ins.GridIndependenceStudy()
#############################################################################
    
########################## Parallel Coordinates #############################  
#    FileName = 'Results.csv'
#    ins.ParallelCoord(FileName)
#############################################################################
    
################################## DOE ######################################     
    N = 900
#    ins.DOE(N)
    lhd = pd.read_csv(ins.Dir+'/Sample.csv').values
    lhd = lhd[423:]
    ins.DOE(N,lhd=lhd)

#############################################################################
    
    

    

    




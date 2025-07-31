# Module import
import os
import sys
import warnings
import time

# Import paths file
pathPath = "/home/s194112/"
sys.path.append(pathPath)
import paths

# Import data folder creation script
sys.path.append(paths.scriptPath)
from create_data_folder import create_data_folder

class bsubmissions:
    """
    A class for submitting jobs into the Sophia HPC dtu uses. All variables need to be filled otherwise a job will fail.
    
    We set jobName, which appears in the console of the jobs
    numNodes -> How many nodes (computers) we wish to use
    tasksPerNode -> How many tasks should be set up per node. Nodes on Sophia consists of 2x16 core processers, so 32 per node are possible
    cpusPerTask -> How many tasks/threads should be run per core. Generally leave this at '1', to ensure single-threading (no wait-time for threads, optimises simulation speed by minimizing wait, 2 threads would optimize use of the cores instead, which we are not interested in)
    scriptName -> What we will call the script we want to submit
    nemPerCPU -> Memory reserved for each core, in GB [Int]. False per default, meaning inactive.
    wallTime -> Max time out program should run for [hh:mm:ss, string], 24:00:00 is default.
    """
    def __init__(self, jobName, scriptName, numNodes = 1, tasksPerNode = 1, cpusPerTask=1, memPerCPU = False, partition = "workq", wallTime = "48:00:00"): # 48 hour timelimit on the workq
        # Set variables
        self.homeFolder = "/home/s194112/" # Set to home directory
        self.processesAsked = int(numNodes*tasksPerNode)
        self.jobName = jobName
        self.nodesAsked = numNodes
        self.tasksPerNode = tasksPerNode
        self.cpusPerTask = cpusPerTask
        self.scriptName = scriptName
        self.partition = partition
        self.wallTime = wallTime
        self.memPerCPU = memPerCPU
    
    # make error paths
    def make_subfolders(self):
        if not (os.path.isdir(self.homeFolder + "/jobs_errors")):
            os.system('mkdir jobs_errors')
        if not (os.path.isdir(self.homeFolder + "/jobs_outputs")):
            os.system('mkdir jobs_outputs')


    def make_sbatch_file(self):
        """ 
        A sbatch file maker, ready to go, fills in all the necessities and creates two subfolders for errors and output files.
        Default behavior is to expect to run the script in the same folder as the submission script, unless an absolute or relative path has been given in the scriptname.
        """
        # Making subfolders for error
        self.make_subfolders() # Are made in home directory

        errorPath = self.homeFolder + '/jobs_errors/'
        outputPath = self.homeFolder + '/jobs_outputs/'
        
        sbatchString = f"""#!/bin/bash
### -- set walltime limit: hh:mm:ss --
#SBATCH --time {self.wallTime}

### -- ask for minimum {self.nodesAsked} number of nodes (can be multiple cpus per node, expect one) (default: 1) --
#SBATCH --nodes {self.nodesAsked}
### -- specify number of tasks/threads per node/pc (2x16 possible per node) (default: 1) --
#SBATCH --ntasks-per-node {self.tasksPerNode}
### -- specify amount of processors per task/thread (keep to 1, then it is one thread per cpu, which is usual for HPCs) (default: 1) --
#SBATCH --cpus-per-task {self.cpusPerTask}

### -- specify partition --
#SBATCH --partition {self.partition}

### -- set the job Name --
#SBATCH --job-name {self.jobName}"""

        # Potentially add max memory per cpu to string
        if self.memPerCPU != False:
            sbatchString = sbatchString + "\n" + f"""
    ### -- specify that we need {self.mem}GB of memory per core/slot --
    #SBATCH --mem-per-cpu {self.mem}"""

        # Keep writing string
        sbatchString = sbatchString + os.linesep + f"""
### -- Specify the output and error file. %J is the job-id --
### -- -o and -e mean append, -oo and -eo mean overwrite -- (is this still true?)
#SBATCH --output {outputPath}Output_{self.jobName}_%J.out 
#SBATCH --error {errorPath}Error_{self.jobName}_%J.err


# here follow the commands you want to execute
cd {self.homeFolder}"""
        if self.processesAsked == 1:
            sbatchString = sbatchString + os.linesep + f"""python3 {self.scriptName}"""
        else:
            sbatchString = sbatchString + os.linesep + f"""mpirun --mca mpi_thread_multiple 0 -np {self.processesAsked} {self.scriptName}
        """
        shFile = open(f'{self.jobName}.sh','w')
        shFile.write(sbatchString) # f"""gammel mpirun -np {self.processesAsked} {self.scriptName} # ny f"""mpirun --mca mpi_thread_multiple 0 -np {self.processesAsked} {self.scriptName}


    def delete_sbatch_file(self):
        os.remove(f'{self.jobName}.sh')

    def do_submission(self, writeJobId = False, dependency = None):
        """
        Possible dependencies are: afterany (starts after it has finished in any state), afterok (after it finishes successfully), etc. See https://bioinformaticsworkbook.org/Appendix/HPC/SLURM/submitting-dependency-jobs-using-slurm.html#gsc.tab=0.
        """
        
        # Construct the string to call
        callStr = ""

        # First, whether the job id must be saved
        if writeJobId == True:
            callStr += f"job_id=$(sbatch --parsable "
        else:
            callStr += f"sbatch "

        # Then, if the job has a dependency
        if dependency == None:
            callStr += f"{self.jobName}.sh"
        else:
            with open(f"{paths.sourcePath}/lastId.txt", 'r') as file:
                prevId = file.read() # read id
            callStr += f"--dependency={dependency}:{prevId} {self.jobName}.sh"

        # Closing parenthesis
        if writeJobId == True:
            callStr += ")"
        
        # Then, closing the string in case we need to save it, and save it (overwriting the file, '>' is for overwriting, '>>' is for appending, in bash)
        if writeJobId == True:
            callStr += f" && echo -n $job_id > ~/lastId.txt"

        # Finally, actually execute the job
        os.system(callStr)
            
        #os.system(f'job_id=$(sbatch --parsable {filename}) ' + rep_no*f'&& job_id=$(sbatch --parsable --dependency=afterok:$job_id {filename}) ')
        #os.system(f'sbatch {self.jobName}.sh')

    def create_and_submit(self, writeJobId = False, dependency = None):
        self.make_sbatch_file()
        self.do_submission(writeJobId = writeJobId, dependency = dependency)
        self.delete_sbatch_file()
        print(f"Jobname '{self.jobName}' has been submitted to the cluster")

    def delete_last_id(self):
        os.remove(f"{paths.scriptPath}/lastId.txt")

class simsubmissions:
    """
    Class for quick submission of nHESEL simulation jobs to the Sophia HPC. Amount of tasks per node is set to 32 (2x16 physical cores per node, for all nodes on Sophia), with 1 'cpu' per task (by cpu they mean processors, so individual cores per task. This here ensure each core runs a single thread per task.)

    inputFile -> Name of input file [str] in HESEL_inputs, used to generate the jobName and submission script
    numNodes -> Amount of nodes/pcs to run on [int]. The total amount of threads to run is nodes * 32, with 2 being the default amount of nodes.
    memPerCPU -> Amount of GB memory to reserve per core [int], left as FALSE by default.
    wallTime -> Max amount of time ["hh:mm:ss"] to run simulation for, before terminating script. Default is "24:00:00".
    overwriteData -> Whether the data folder being created should be overwritten, if already present. True/False, default: False.
    """

    def __init__(self, inputFile, numNodes = 2, memPerCPU = False, wallTime = "24:00:00", overwriteData = False):
        self.inputFile = inputFile
        self.numNodes = numNodes
        self.memPerCPU = memPerCPU
        self.wallTime = wallTime
        self.overwriteData = overwriteData

        # Generate path for data file & source path
        self.dataFolderPath = paths.dataPath + f"/data_{inputFile}/" # Doing it here allows for overwriting self.dataFolderPath with some other path, and then running the rest of the script
        self.HESELPath = paths.sourcePath + "/BOUT-HESEL/hesel"

    def create_data_folder(self):
        # Create a data folder if not present, or overwrite a present folder if overwriteData is set to True
        create_data_folder(inputName = self.inputFile, overwrite = self.overwriteData)
    
    def generate_jobName(self):
        # Generate jobname from input file-name
        self.jobName = "nHESEL_" + self.inputFile

    def generate_scriptPath(self):
        # Generate hesel simulation submission command
        self.HESELsub = self.HESELPath + " -q -d " + self.dataFolderPath

    def do_submission(self, writeJobId = False, dependency = None):     
        # Generate bsubmission object for job submission
        submission = bsubmissions(jobName = self.jobName, scriptName = self.HESELsub, numNodes = self.numNodes, tasksPerNode = 32, memPerCPU = self.memPerCPU, wallTime = self.wallTime)
        
        # Do submission
        submission.create_and_submit(writeJobId = writeJobId, dependency = dependency)

    def do_post_sim_submission(self, jobName, scriptName, writeJobId = False, dependency = None):
        # Generate bsubmission object for job submission
        submission = bsubmissions(jobName = jobName, scriptName = scriptName, numNodes = 1, tasksPerNode = 1, memPerCPU = self.memPerCPU, partition = "workq", wallTime = "24:00:00")
        
        # Do submission
        submission.create_and_submit(writeJobId = writeJobId, dependency = dependency)

    def ready_and_submit(self):
        # Do the whole usual routine for submitting a sim job
        # First ensure the data folder is present, and if so, give a warning if overwriteData is set to false and continue without creating a new folder
        if os.path.isdir(self.dataFolderPath) == True and self.overwriteData == False:
            warnings.warn(f"The indicated input file '{self.inputFile}' already has a data-folder present at '{self.dataFolderPath}', and overwriteData is set to '{self.overwriteData}'. Folder already present will be used in simulation job.")
        # Otherwise simply create or overwrite the data folder
        else:
            self.create_data_folder()

        # Now generate the jobname and scriptpath
        self.generate_jobName()
        self.generate_scriptPath()

        # Finally do the submission of the simulation job
        self.do_submission()

    def do_full_process(self, includeSim = True, includeProcessing = True):
        """
        Usual routine for submitting a sim job, the postprocessing, and default figures.

        Note that the data folders MUST be on the format data_{input file name}.

        The simulation and processing steps can be turned off by setting includeSim and includeProcessing, respectively, to False.

        """

        if includeSim == True:
            # First ensure the data folder is present, and if so, give a warning if overwriteData is set to false and continue without creating a new folder
            if os.path.isdir(self.dataFolderPath) == True and self.overwriteData == False:
                warnings.warn(f"The indicated input file '{self.inputFile}' already has a data-folder present at '{self.dataFolderPath}', and overwriteData is set to '{self.overwriteData}'. Folder already present will be used in simulation job.")
            # Otherwise simply create or overwrite the data folder
            else:
                self.create_data_folder()
        else:
            if os.path.isdir(self.dataFolderPath) != True:
                raise Exception(f"No datafolder is present at {self.dataFolderPath} for processing or figures. Aborting.")

        if includeSim == True:
            ## First the sim submission
            # Now generate the jobname and scriptpath for the sim
            self.generate_jobName()
            self.generate_scriptPath()
    
            # Submit the sim, saving the id
            self.do_submission(writeJobId = True, dependency = None)
    
            # 'lil wait
            time.sleep(1)

        ## Now the post-processing submissions
        ## First the processing
        if includeProcessing == True:
            # Start by writing using a template file
            processingFileName = f"{paths.tempPath}/processingTemp_{self.inputFile}.py"
            with open(f"{paths.templatePath}/TEMPLATE_processing.py","r") as templateFile:
                tempFile = open(processingFileName,"w")
                for line in templateFile:
                    updatedLine = line.replace("NAMEDATA", f"data_{self.inputFile}")
                    tempFile.write(updatedLine)
                tempFile.close()
    
            # Then submit it
            if includeSim == True:
                self.do_post_sim_submission(jobName = f"{self.inputFile}-Processing", scriptName = processingFileName, writeJobId = True, dependency = "afterany")
            else:
                self.do_post_sim_submission(jobName = f"{self.inputFile}-Processing", scriptName = processingFileName, writeJobId = True)
    
            # 'lil wait
            time.sleep(1)
        else:
            pass

        # Then delete it again
        #os.remove(f"{paths.sourcePath}/processingTemp.py") # This needs to be done as a seperate job *after* the processing job is done, as it otherwise results in the job being unable to execute, due to the file *not existing at execution*

        ## Then the figures
        # Start by writing using a template file
        figuresFileName = f"{paths.tempPath}/figuresTemp_{self.inputFile}.py"
        with open(f"{paths.templatePath}/TEMPLATE_plotting.py","r") as templateFile:
            tempFile = open(figuresFileName,"w")
            for line in templateFile:
                updatedLine = line.replace("NAMEDATA", f"data_{self.inputFile}")
                tempFile.write(updatedLine)
            tempFile.close()

        # Then submit it
        if includeProcessing == True:
            self.do_post_sim_submission(jobName = f"{self.inputFile}-Figures", scriptName = figuresFileName, writeJobId = False, dependency = "afterok")
        else:
            self.do_post_sim_submission(jobName = f"{self.inputFile}-Figures", scriptName = figuresFileName, writeJobId = False)

        # 'lil wait
        time.sleep(1)

        # Then delete it again
        #os.remove(f"{paths.sourcePath}/figuresTemp.py") # Same issue as above

        


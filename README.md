# Matlab-Whisker-Tracker
Welcome to the Matlab Whisker Tracker repository! This toolbox contains functions to extract, analyse and visualize whiskers in video recordings. To get started follow the [quickstart](#quickstart) below. The toolbox documentation is available at the [wiki page](https://github.com/DepartmentofNeurophysiology/Matlab-Whisker-Tracker/wiki).

The MWT was build and tested in MATLAB 2016a and 2017a and requires the 'image processing' and 'signal processing' toolbox.

<p align="center">
   <img  src="https://github.com/DepartmentofNeurophysiology/Matlab-Whisker-Tracker/wiki/Gifs/README.gif"
      width=500/> 
</p> 
 
 
 # Quickstart
 ### Install
 Clone this repository by downloading and extracting the repo's zips package or use git install:
 
 `git clone htps://github.com/DepertmentofNeurophysiology/Matlab-Whisker-Tracker  LOCAL_DIRECTORY`
 
 Add the repository to the matlab path:
 
 `addpath(genpath(LOCAL_DIRECTORY))`
 
 ### First time setup
 Verify the tracker settings by reading trough './Settings/defaultSettings.m' and update variables where needed. When finished save and close the file.
 
 Run `Parametersetup` to generate a './Settings/Settings.mat' file used for video tracking. The functions allows the user to finetune tracker parameters.
 Select a folder containing videos for tracking (it searches subdirectories as well), when a folder is selected then panel depicted below shows up. On the left side information of the data displayed is visible: Current path (data path selected), Video Select (a list of found videos, and the video currently on display) and Current frame (the current frame on display). Change the displayed video or selected frame using the buttons.
 Each panel shows a different processing step and has sliders to interact with relevant parameters. When hovering over the slidernames, tooltips describe the slider target variable:
 * GK.      - Gaussian kernel size
 * BG. t.   - Background threshold
 * EG. ks.  - Delta kernel size for edge detection
 * EG. kl.  - Box kernel size for edge detection
 * EG. t.   -  Edge threshold
 * SE. t.   -  Shape threshold
 * SE. dil. - Dilation size of ROI
 * SD. t.   - Seed threshold
 * TT. t.   - Trace threshold
 * TT. ks.  - Delta kernel size for trace detection
 * TT. kl.  - Box kernel size for edge detection
 
 Further check boxes can be used to change displayed data. Use the save button to save data, and the quit button to exit the GUI. The presented results in the ParameterSetup are the same output as from the full tracking algorithm.
 
<p align="center">
<img  src="https://github.com/DepartmentofNeurophysiology/Matlab-Whisker-Tracker/wiki/Images/ParameterSetup.png"
      width=550/>
</p>
 
 ### Video tracking
 Run Video tracking by executing `BatchTracking`. A matlab dialog will ask for a folder containing videofiles. After selecting a a video folder, a GUI shows up with all detected videos. The column on the left, 'Videos to track', contains all videos that will be tracked after pressing 'confirm'. The column on the right, 'Videos to skip', contains all videos that will not be tracked. 
 
 Select videos in 'Videos to track' and press 'Remove' to remove the videos from the job list. Select videos in 'Videos to skip' and press 'Add' to add videos to the job list. By default, videos that have not been tracked before are listed under 'Videos to track' and videos that have been tracked before are listed under 'Videos to skip' and marked with an exclamation-mark (!). **Press confirm** to initiate tracking.
 <p align="center">
   <img src="https://github.com/DepartmentofNeurophysiology/Matlab-Whisker-Tracker/wiki/Images/BatS1.png"
        width =500/>
 </p>
 
 ### Output
 The tracker data is saved in the video path with the name 'VIDEO_NAME_Annotations_Tracker.mat'. The returned output is described in the [Data-storage](https://github.com/DepartmentofNeurophysiology/Matlab-Whisker-Tracker/wiki/Data-storage) page. For reading on further analysis and visualization, read the repository [wiki](https://github.com/DepartmentofNeurophysiology/Matlab-Whisker-Tracker/wiki).
 
 


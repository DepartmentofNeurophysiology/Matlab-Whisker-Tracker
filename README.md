# Matlab-Whisker-Tracker
Welcome to the Matlab Whisker Tracker repository! This toolbox contains functions to extract, analyse and visualize whiskers in video recordings. To get started follow the [quickstart](#quickstart) below. The toolbox documentation is available at the [wiki page](https://github.com/DepartmentofNeurophysiology/Matlab-Whisker-Tracker/wiki).

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
 Verify the tracker settings by reading trough './Settings/makeSettings.m' and update variables where needed. When finished close and save the file.
 
 Run `ParameterSetup`, the function will show a GUI to select a threshold for background extraction. **Update the current threshold** by moving the sliderbar. When a proper threshold is found, **press 'confirm'**. Next, a GUI is shown to select a frame for demo-tracking. Use the sliderbar to **select a frame** to be tested and press confirm. Now a final GUI is shown with tracking results for the selected frame, with the current tracker settings. Three parameters can be updated manually:
 
 * Dilation size - distance of tracing seeds from the head
 * Origin threshold (OR tr.) - sensitivity of seed detection
 * Tracking threshold (TR tr.) - threshold to stop tracking a trace
 
 Press 'overlay' to visualize the edges of detected objects, 'switch frame' to select another frame. After settings the correct parameters, press **'confirm'** to save the updated settings and quit the tracker setup.
<p align="center">
<img  src="https://github.com/DepartmentofNeurophysiology/Matlab-Whisker-Tracker/wiki/Images/ParS1.png"
      width=250/>
<img  src="https://github.com/DepartmentofNeurophysiology/Matlab-Whisker-Tracker/wiki/Images/ParS2.png"
      width=250/>
<img  src="https://github.com/DepartmentofNeurophysiology/Matlab-Whisker-Tracker/wiki/Images/ParS3.png"
      width=250/>
</p>
 
 ### Video tracking
 Run Video tracking by executing `BatchTracking`. A matlab dialog will ask for a folder containing videofiles. After **selecting a a video folder**, a GUI shows up with all detected videos. The column on the left, 'Videos to track', contains all videos that will be tracked after pressing 'confirm'. The column on the right, 'Videos to skip', contains all videos that will not be tracked. 
 
 Select videos in 'Videos to track' and press 'Remove' to remove the videos from the job list. Select videos in 'Videos to skip' and press 'Add' to add videos to the job list. By default, videos that have not been tracked before are listed under 'Videos to track' and videos that have been tracked before are listed under 'Videos to skip' and marked with an exclamation-mark (!). **Press confirm** to initiate tracking.
 <p align="center">
   <img src="https://github.com/DepartmentofNeurophysiology/Matlab-Whisker-Tracker/wiki/Images/BatS1.png"
        width =500/>
 </p>
 
 ### Output
 The tracker data is saved in the video path with the name 'VIDEO_NAME_Annotations_Tracker.mat'. The returned output is described in the [Data-storage](https://github.com/DepartmentofNeurophysiology/Matlab-Whisker-Tracker/wiki/Data-storage) page. For reading on further analysis and visualization, read the repository [wiki](https://github.com/DepartmentofNeurophysiology/Matlab-Whisker-Tracker/wiki).
 
 


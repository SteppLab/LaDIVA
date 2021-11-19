# LaDIVA
## DIVA connected to BCM Zañartu 2014
### V007 - Hardcorded for Pitch Reflex : run one iteration, change perturbation magnitude, duration, direction
        * BCMlookuptable.m download from https://drive.google.com/file/d/1Arn8HqNxDIot3IpTC7q68Wmqj7G_--lC/view?usp=sharing
        * Current version upladed to PLOS Comp Bio submission
        
        # LaDIVA sourcecode
The source code of Laryngeal DIVA model, which combines DIVA model of speech motor contorl with extended body cover model for laryngeal biomechanical control. 

## Installation
1. Download the sourcecode provided to a target location (e.g.: localpath/LaDIVA) in the local machine.
2. Download the mat file 'BCMlookuptable.m' from the link [here] and add it to the same target folder (localpath/LaDIVA/) the sourcecode was saved to. 

[here]: https://drive.google.com/file/d/1Arn8HqNxDIot3IpTC7q68Wmqj7G_--lC/view

## Running the simulink model
1. Click on the model file 'diva.slx' under the target folder (localpath/LaDIVA/diva.slx) to open the simulink model and click on the **Production 'input'** module to open the DIVA Graphic User Interface(GUI).
![Fig_divaSimulink](https://user-images.githubusercontent.com/13642912/138382798-71486928-02d1-4485-9537-185726db9bcf.PNG)
<p align = "left"> <b>Figure 1.LaDIVA Simulink Model</b></p>

2. Select the **input** and **output** signals that needs to be observed. In the example below Cricothyroid(CT) muscle activation, Thyroarytenoid(TA) muscle activation, Subglottal pressure, and glottal constriction signals are plotted as input signals, and vocal fundermental frequency (*f*o), vocal Sound Pressure Level (SPL), and voicing are plotted as output signals. 

![Fig_DivaGUI](https://user-images.githubusercontent.com/13642912/138383548-401f7220-cdfb-4b45-893d-84cbc8fbb25c.PNG)
<p align = "left"> <b>Figure 2.LaDIVA Graphic User Interface(GUI)</b></p>

3. Select the target sound that needs to be produced from the drop down menu at the bottom and click **Start** button to run a simulation.

## Introducing perturbations to the model
1. Auditory reflexive or auditory adaptive perturbations of vocal *f*o can be introduced by navigating to the simulink model and opening the **auditory perturbation** block.
2. First select the Perturbation type
3. For Auditory Reflexive Perturbation, you can run one reflexive trial simulation. Set perturbation magnitude (in Cents), perturbation duraiton, and perturbation onset as propted. See Fig 3.

![PitchReflex](https://user-images.githubusercontent.com/13642912/142690571-8db28e69-55a9-40dd-9f47-7879211a7dfa.JPG)
<p align = "left"> <b>Figure 3. For Pitch Reflexive Perturbation: Setting perturbation magnitude, perturbation duration, and perturbation onset in Auditory Perturbation block</b></p>

4. For Auditory Adaptive Perturbation, you can run multiple trial length simulations. Set maximum perturbation magnitude (in Cents), and number of trials in each phase as propted. See Fig 4. 
*Note: you should select the number of trials as the number of repeatitions in DIVA_GUI*

![PitchAdapt](https://user-images.githubusercontent.com/13642912/142690702-048cd25a-51f2-4f52-9b17-56e3c1e6cb76.JPG)
<p align = "left"> <b>Figure 4. For Pitch Adaptive Perturbation: Setting Phases in adaptation paradigm and max perturbation magnitude in Auditory Perturbation block</b></p>

## Modifying Feedback and Feedforward Gains
1. In the main Simulink model > select 'Articulatory Velocity and Position Map" > and click on Gain blocks to modify feedback and feedforward gains

![GainControl](https://user-images.githubusercontent.com/13642912/142690174-85eea060-2c5b-49da-8834-c2fb3088da33.JPG)
<p align = "left"> <b>Figure 5. Artucilatory Velocity and Position Map block</b></p>

3. You can modify the following gain values
   * Feedback Gain (overall) : This is set to 1 and not modified
   * Feedforwward Gain (overall) : This is set to 1 and not modified
   * Auditory Feedback Gain : Modify to fit auditory reflexive experiment datasets
   ![FBgain](https://user-images.githubusercontent.com/13642912/142690287-8a008cd6-c4ac-4512-b8f3-330b64bca77e.JPG)
   
   * Somatosensory Feedback Gain : This is set to 0 and not modified                
   ![FBgainsom](https://user-images.githubusercontent.com/13642912/142690309-c609400c-872a-4bd0-94d3-9f197908b06f.JPG)
   
   * Feedforward Learning Rate : Modify to fit auditory adaptive experiment datasets
   ![learningrate](https://user-images.githubusercontent.com/13642912/142690324-ab3050e2-80fb-4d36-ab61-af41462dcf98.JPG)


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
<p align = "center"> <b>Figure 1.LaDIVA Simulink Model</b></p>

2. Select the **input** and **output** signals that needs to be observed. In the example below Cricothyroid(CT) muscle activation, Thyroarytenoid(TA) muscle activation, Subglottal pressure, and glottal constriction signals are plotted as input signals, and vocal fundermental frequency (*f*o), vocal Sound Pressure Level (SPL), and voicing are plotted as output signals. 

![Fig_DivaGUI](https://user-images.githubusercontent.com/13642912/138383548-401f7220-cdfb-4b45-893d-84cbc8fbb25c.PNG)
<p align = "center"> <b>Figure 2.LaDIVA Graphic User Interface(GUI)</b></p>

3. Select the target sound that needs to be produced from the drop down menu at the bottom and click **Start** button to run a simulation.

## Introducing perturbations to the model
1. Auditory reflexive or auditory adaptive perturbations of vocal *f*o can be introduced by navigating to the simulink model and opening the **auditory perturbation** block.

![Fig_perturbationType](https://user-images.githubusercontent.com/13642912/138384876-90eb1312-9236-464a-9a55-15e74458f8b6.PNG)
<p align = "center"> <b>Figure 3.Selecting perturbation type, magnitude, and duration from Auditory Perturbation block</b></p>

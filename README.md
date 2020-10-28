# SEIR-DIT
SEIR-DIT is a mathematical compartmental model aimed to study the effect of DIT strategy (Detect symptoms, isolate and trace contacts) as an alternative to COVID-19 control.
The model is composed of two main ODE systems. A SEIR-like model (see Figure 1) that describes the epidemic dynamics and contains the main epidemiological states, and a secondary compartmental model (see Figures 2,3 and 4), fed by the main SEIR-like model, that allows to calculate the testing requirements and amount of isolated cases and contacts at every time, necessary for the implementation of the strategy. 

![Esquema](/IMG/SEIR-MODEL.png?raw=true)
##### Fig. 1. SEIR-like model
<p></p>
<p></p>
![Esquema](/IMG/Quarantine-Probable.png?raw=true)
##### Fig. 2. Quarantine and testing of probable index cases 
<p></p>
<p></p>
![Esquema](/IMG/Quarantine-SuspectedAndContacts.png?raw=true)
##### Fig. 3. Quarantine and testing of probable index cases and their contacts
<p></p>
<p></p>
![Esquema](/IMG/Quarantine-Exposed.png?raw=true)
##### Fig. 4. Quarantine and testing of exposed contacts 
<p></p>
<p></p>
![Esquema](/IMG/Quarantine-NonExposed.png?raw=true)
##### Fig. 3. quarantine and testing of non-exposed contacts 







# SEIR-DIT
SEIR-DIT is a mathematical compartmental model aimed to study the effect of DIT strategy (Detect symptoms, isolate and trace contacts) as an alternative to COVID-19 control.
The model is composed of two main ODE systems. A SEIR-like model (see Figure 1) that describes the epidemic dynamics and contains the main epidemiological states, and a secondary compartmental model (see Figures 2,3 and 4), fed by the main SEIR-like model, that allows to calculate the testing requirements and amount of isolated cases and contacts at every time, necessary for the implementation of the strategy. 

## Main SEIR-like model
<ul>
  <li>S: Susceptible</li>
  <li>Q<sub>S1</sub>: Quarantine for suspected index cases and contacts</li>
  <li>Q<sub>S2</sub>: Quarantine for non-exposed traced contacts of detected mild infections</li>
  <li>E: Exposed</li>
  <li>E<sup>T</sup>: Exposed traced contacts of  detected mild infections</li>
  <li>I<sub>M</sub>: Mild infections</li>
  <li>I<sub>M</sub><sup>D</sup>: Detected mild infections</li>
  <li>I<sub>M</sub><sup>T</sup>: Traced mild infections (coming from exposed traced contacts)</li>
  <li>I<sub>C</sub>: Infections requiring hospitalization</li>
  <li>I<sub>C</sub><sup>SI</sup>: Self-isolated infections requiring hospitalization</li>
  <li>I<sub>C</sub><sup>T</sup>: Traced infections requiring hospitalization (coming from exposed traced contacts) </li>
  <li>I<sub>HR</sub>: Hospitalized infections requiring a general hospital bed that recover</li>
  <li>I<sub>UR</sub>: Hospitalized infections requiring an ICU bed that recover</li>
  <li>I<sub>HD</sub>: Hospitalized infections requiring a general hospital bed that die</li>
  <li>I<sub>UD</sub>: Hospitalized infections requiring an ICU bed that die</li>
  <li>I<sub>R</sub>: Hospitalized infections requiring a general hospital bed after recovering from ICU stay</li>
  <li>R: Recovered</li>
  <li>D: Dead</li>
</ul>


|![Esquema](/IMG/SEIR-MODEL.png?raw=true)|
| :---:         | 
|**Fig. 1.** SEIR-like model|

## Secondary compartmental model for testing and isolation
<ul>
  <li>Q<sub>IM</sub><sup>D1</sup>: Detected mild infections with positive antigen test</li>
  <li>Q<sub>S1</sub>: Quarantine for suspected index cases and contacts</li>
  <li>Q<sub>S2</sub>: Quarantine for non-exposed traced contacts of detected mild infections</li>
  <li>Q<sub>E</sub><sup>T</sup>: Exposed traced contacts (isolated after tracing)</li>
  <li>Q<sub>IM</sub><sup>T</sup>: Isolated mild infections (coming from traced exposed contacts)</li>
  <li>Q<sub>IM</sub><sup>T1</sup>: Traced mild infections with positive RT-PCR test</li>
  <li>Q<sub>IC</sub><sup>T</sup>: Isolated infections requiring hospitalization (coming from traced exposed contacts)</li>
</ul>

| ![Esquema](/IMG/Quarantine-Probable.png?raw=true) | ![Esquema](/IMG/Quarantine-SuspectedAndContacts.png?raw=true) |
| :---:         |     :---:      | 
| **Fig. 2.** Isolation and testing of probable index cases | **Fig. 3.** Isolation and testing of suspected index cases and their contacts |

| ![Esquema](/IMG/Quarantine-Exposed.png?raw=true) |
| :---:         | 
| **Fig. 4.** Isolation and testing of exposed contacts |

|![Esquema](/IMG/Quarantine-NonExposed.png?raw=true)|
| :---:         |
|**Fig. 5.** Isolation and testing of non-exposed contacts |

## Parameters
### Epidemiological parameters for transitions between infectious states
The parameters of the main SEIR-like model were adapted from https://mrc-ide.github.io/global-lmic-reports/parameters.html 
|Parameter | Value | Definition |
|---|---|---|
|R0|3.0|Reproduction number|
|β<sub>0</sub>|1.3743 |Transmission rate of infectious individuals that move freely |
|β<sub>1</sub>|0.0199 |Transmission rate of infectious individuals in households |
|b| - |Probability of transmitting the disease b=β<sup>M</sup>/n |
|n|15.721|Average number of contacts per individual|
|β<sup>H</sup>|0.001|Hospital transmission rate|
|ω<sup>-1</sup> |4.6 days|Mean incubation period |
|σ<sub>C</sub><sup>-1</sup>|3 days|Mean time upon self-isolation for severe infections|
|σ<sub>CSI</sub><sup>-1</sup>|4.1 days|Mean duration of isolation for severe infections prior hospitalization|
|σ<sub>CT</sub><sup>-1</sup>|7.1 days|Mean duration of isolation for traced severe infections prior hospitalization|
|γ<sub>M</sub><sup>-1</sup>|2.1 days|Mean duration of mild infection|
|γ<sub>HR</sub><sup>-1</sup>|9.5 days|Mean duration of hospitalization for non-critical cases if survive|
|ν<sup>-1</sup>|11.3 days|Mean duration in ICU if survive|
|γ<sub>R</sub><sup>-1</sup>|3.4 days|Mean duration of stepdown post ICU|
|σ<sub>HD</sub><sup>-1</sup>|7.6 days|Mean duration of hospitalization for non-critical cases if die|
|σ<sub>UD</sub><sup>-1</sup>|10.1 days|Mean duration in ICU if die|
|δ<sub>M</sub>|0.96558|Probability of mild infections |
|δ<sub>HR</sub>|0.69660|Probability of hospitalized infections requiring a general hospital bed that recover from disease|
|δ<sub>UR</sub>|0.12257|Probability of hospitalized infections requiring an ICU bed that recover |
|δ<sub>HD</sub>|0.05827|Probability of hospitalized infections requiring a general hospital bed that die|
|δ<sub>UD</sub>|0.12256|Probability of hospitalized infections requiring an ICU bed that die|

### DIT parameters 
#### A. Contact tracing

|Group|Contacts per group|SAR|Proportion of traced contacts |Traced contacts of suspected index cases (n<sub>T</sub>)|\[\phi\]|a|
|---|---|---|---|---|---|---|

Index case
0
-
-
0
0
0
Household
3.039
25%
100%
3.039
0.37
0.12
Work / school
5.931
15%
80%
7.784
0.74
0.42
Other
6.751
5%
50%
11.160
0.81
0.67





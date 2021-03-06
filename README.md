# SEIR-DIT
SEIR-DIT is a mathematical compartmental model aimed to study the effect of detection of symptoms, isolation and contacts tracing as an alternative to COVID-19 control.

The model is composed of two main ODE systems. A SEIR-like model (see Figure 1) that describes the epidemic dynamics and contains the main epidemiological states, and a secondary compartmental model (see Figures 2,3 and 4), fed by the main SEIR-like model, that allows to calculate the testing requirements and the amount of isolated cases and contacts at every time. 

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
  <li>I<sub>HR</sub>: Hospitalized infections that require a general hospital bed and recover</li>
  <li>I<sub>UR</sub>: Hospitalized infections that require an ICU bed and recover</li>
  <li>I<sub>HD</sub>: Hospitalized infections that require a general hospital bed and die</li>
  <li>I<sub>UD</sub>: Hospitalized infections that requiring an ICU bed and die</li>
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
The parameters of the main SEIR-like model were adapted from https://mrc-ide.github.io/global-lmic-reports/parameters.html, assuming that the probability of death for critical cases is 0.35
|Parameter | Value | Definition |
|---|---|---|
|R0|3.0|Reproduction number|
|β<sub>0</sub>|1.3743 |Transmission rate in absence of social distancing |
|β<sub>1</sub>|0.0199 |Transmission rate within households |
|b| - |Probability of transmitting the disease b=β<sub>M</sub>/n |
|n|15.721|Average number of contacts per individual|
|β<sub>H</sub>|0.001|Hospital transmission rate|
|ω<sup>-1</sup> |4.6 days|Mean incubation period |
|σ<sub>C</sub><sup>-1</sup>|3 days|Mean time upon self-isolation for severe infections|
|σ<sub>CSI</sub><sup>-1</sup>|4.1 days|Mean duration of isolation for severe infections prior hospitalization|
|σ<sub>CT</sub><sup>-1</sup>|7.1 days|Mean duration of isolation for traced contacts with severe infection prior hospitalization|
|γ<sub>M</sub><sup>-1</sup>|2.1 days|Mean duration of mild infection|
|γ<sub>HR</sub><sup>-1</sup>|9 days|Mean duration of hospitalization for non-critical cases if survive|
|ν<sup>-1</sup>|14.8 days|Mean duration in ICU if survive|
|γ<sub>R</sub><sup>-1</sup>|3 days|Mean duration of stepdown post ICU|
|σ<sub>HD</sub><sup>-1</sup>|9 days|Mean duration of hospitalization for non-critical cases if die|
|σ<sub>UD</sub><sup>-1</sup>|11.1 days|Mean duration in ICU if die|
|δ<sub>M</sub>|0.97902 |Probability of mild infection|
|δ<sub>HR</sub>|0.67959|Probability of recovery for hospitalized infections requiring a general hospital bed|
|δ<sub>UR</sub>|0.13883|Probability of recovery for hospitalized infections requiring an ICU bed|
|δ<sub>HD</sub>|0.10682|Probability of dying for hospitalized infections requiring a general hospital bed|
|δ<sub>UD</sub>|0.07476|Probability of dying for hospitalized infections requiring an ICU bed |

### DIT parameters 
#### A. Contact tracing

  |Group|Contacts per group|SAR|Proportion of traced contacts |Traced contacts of suspected index cases (n<sub>T</sub>)|![formula](https://render.githubusercontent.com/render/math?math=\Phi)|![formula](https://render.githubusercontent.com/render/math?math=\tilde\Phi)|
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
|Index case|0|-|-|0|0|0|
|Household|2.039|25%|100%|2.039|0.26|0.11|
|Work / school|5.931|15%|80%|6.784|0.61|0.41|
|Other|7.751|5%|50%|10.659|0.71|0.67|

#### B. Testing and isolation
|Parameter|Value|Definition|
|:---:|:---:|:---|
|α<sup>-1</sup>|1 day|Mean time between onset of symptoms and detection|
|ρ<sup>-1</sup>|12 days|Mean duration of isolation period|
|T<sub>PCR</sub>|2 days|Mean time to results of RT-PCR test|
|ξ<sub>PCR</sub>|0.85|Sensitivity of RT-PCR test|
|T<sub>AG</sub>|1 days|Mean time to results of AG test|
|ξ<sub>AG</sub>|0.75|Sensitivity of AG test|

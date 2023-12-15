![image](https://github.com/MiltonCardenas/Grad_Math2/assets/121530690/23256736-f539-4350-8a81-dc891ba8bf2f)

# Modeling and Estimation of Kinetic Parameters and Replicative Fitness of HIV-1
The objective of the present GitHub Repository is to analyze how HIV-1 replicates by using the dynamic models proposed in [1]. Especially, this analyzis involves 3 main components: data fitting of experimental information, sensitivity analysis and a simplified bifurcation analysis. 

Architecture:
The majority of the codes presented are Python files, created in the Spyder scientific environment. Results are summarized in graphics within the same files as the original code.

# Introduction
HIV-1, or Human Immunodeficiency Virus Type 1, is the primary cause of HIV/AIDS in humans, a virus that targets the immune system, particularly CD4+ T cells, leading to their gradual depletion. This weakening of the immune system leaves individuals vulnerable to a wide range of opportunistic infections and diseases. HIV-1 is characterized by its high mutation rate, which results in a diverse population of viral variants within infected individuals. 

## Importance of understanding HIV-1 replication ##
The modeling of HIV-1's replicative dynamics is importnt because it offers a deeper understanding of how the virus adapts and survives in various environments, particularly in response to antiviral treatments. By accurately simulating these dynamics, is possible to predict the behavior of different viral strains, enabling the identification of the most effective treatment strategies. This modeling is essential for tailoring patient-specific therapies and enhancing the overall efficacy of HIV treatment regimens. Furthermore, it aids in anticipating the virus's evolution, potentially leading to the development of more effective drugs and vaccines. In essence, these models are a vital tool in the ongoing battle against HIV, providing insights that can lead to improved patient outcomes and a better understanding of the virus's complex behavior.

# Objectives

- Estimate the necessary parameters for modeling the HIV-1 Replicative Fitness.
- Do a Sensitivity analysis of the Model by modifying some key parameters.
- Find the steady states of the HIV-1 replicative dynamics.

# Dynamic Model
The following equations describe the dynamics of cell growth and virus infection, this model was taken from [1]:

- Cells growth equation:
  $$\frac{dT}{dt} = (\rho - k_m T_m - k_w T_w - k_R T_{mw})T$$

- Cells infected by mutant Virus:
  $$\frac{dT_m}{dt} = (\rho_m + k_m T - q_m T_w) T_m + 0.25k_R T_{mw}T$$

- Cells infected by wild-type Virus:
  $$\frac{dT_w}{dt} = (\rho_w + k_w T - q_w T_m) T_w + 0.25k_R T_{mw}T$$

- Cells infected by both types of virus:
  $$\frac{dT_{mw}}{dt} = (\rho_{mw} + 0.5k_R T) T_{mw} + (q_m + q_w)T_w T_m$$

# Dataset

The dataset for the current study consists of three experimental runs focused on growth measurement. Given the inherent challenges associated with measuring these variables (Number of Cells), the collected data was analyzed in conjunction with various growth models to ensure robust correlation and validation.

![image](https://github.com/MiltonCardenas/Grad_Math2/assets/121530690/dfa583ec-f1e3-4a5b-9918-6478ff6247a8)

# Model Fitting

The model fitting consists of 4 ODEs and 9 parameters:

$rho$ $rho_{m}$  $rho_{w}$  $rho_{mw}$  are cell growth rates for healthy, mutant infected, wild infected, and double infected respectively.

$k_{m}$  $k_{w}$ $k_{R}$ are the infection rates for mutant strain, wild strain and the combined infection respectively.

$q_{m}$ $q_{w}$ are dual infection rates for mutant strain and wild strain respectively

To fit the model it was considered:
1 - Test the ODE function with the parameters of the cited paper
2 - Evaluation of the 'minimize' function's performance by assessing its error rates and computational time efficiency.
3 - Use of a logarithmic scale to achieve system convergence and to effectively manage the large disparities among dependent variables in their initial states.
4. Analyze the impact of initial guesses to fitting outcomes.
5. Compare the results

## Model Fitting Results:

Dynamic Profiles of Cells growing (Current work)
![Fitted dynamic model](https://github.com/MiltonCardenas/Grad_Math2/assets/121530690/aebd1c3d-acfb-40cf-8855-c539ae4549f4)

Dynamic Profiles of cells growth (From the results of the paper)
![image](https://github.com/MiltonCardenas/Grad_Math2/assets/121530690/d3673384-3e5e-4df2-a098-c9a5b28c47cd)

Qualitatively, it can be appreciated that the models exhibit similar profiles among themselves, which is a positive result as it reinforces confidence in the set of parameters obtained.

Numeric results are presented in the next image:
![image](https://github.com/MiltonCardenas/Grad_Math2/assets/121530690/90db6813-647c-4e9f-bc7e-17f96518bc32)

To ascertain whether the system's outcomes are influenced by the initial parameter estimates, the initial guesses were varied within a defined threshold from the values reported in the literature. The effects of these variations on the results are summarized in the tables below:

![image](https://github.com/MiltonCardenas/Grad_Math2/assets/121530690/50984288-c85f-4a8f-9651-1819c55f13fa)

It is important to verify the initial parameter estimates through an in-depth analysis. An increase in the number of parameters tends to complicate the model and can detract from the reliability of the results. Therefore, in instances where there is a considerable variance among the parameters, it is advisable to implement normalization techniques or logarithmic transformations.

# Model Sensitivity Analysis

Parameter 1: $\rho$ (Healthy cells growth rate) 
![Normalized sensitivity analysis Param 2](https://github.com/MiltonCardenas/Grad_Math2/assets/121530690/961fbada-f80d-4816-86a1-3697ace274ce)

Parameter 2: $\rho_{m}$ (Mutant strain infected cells growth rate) 
![Normalized sensitivity analysis Param 2](https://github.com/MiltonCardenas/Grad_Math2/assets/121530690/269021eb-551b-4600-8ad9-e7d3df2c6e5c)

Parameter 3: $\rho_{w}$ (Wild strain infected cells growth rate) 
![Normalized sensitivity analysis Param 3](https://github.com/MiltonCardenas/Grad_Math2/assets/121530690/74e395d4-aab2-47ea-98d3-704ba940ec5a)

Parameter 4: $\rho_{mw}$ (Dual infected cells growth rate) 
![Normalized sensitivity analysis Param 4](https://github.com/MiltonCardenas/Grad_Math2/assets/121530690/0de6ff25-541a-4b09-9db2-dba0b7cd8cb4)

The growth rate of healthy cells is the parameter with the most significant impact on the system. Alterations to this parameter substantially influence the sensitivity of the model, particularly concerning the cells infected by the mutant strain, as demonstrated by the next sensitivity analysis:

![Normalized sensitivity analysis Param 1](https://github.com/MiltonCardenas/Grad_Math2/assets/121530690/ef0aa39c-8388-446d-b2ce-f4c11ad798bd)

Sensitivity analysis showcases the impact of parameter variability on the outcomes of the model, underscoring the essential factors that affect the behavior of the system. This analysis is particularly beneficial as it highlights the most critical parameters that should be targeted for data collection. By focusing on these key parameters, researchers can allocate resources more efficiently and improve the accuracy of the model in simulating real-world scenarios.

# Simplified Model Bifurcation Analysis

For the bifurcation analysis, the model was simplified:

- Cells growth equation:
  $$\frac{dT}{dt} = (\rho - k_m T_m) \cdot T$$

- Cells infected by mutant Virus:
  $$\frac{dT_m}{dt} = (\rho_m + k_m T) \cdot T_m$$

Steady states were derived:

$$0 = (\rho - k_m T_m) \cdot T$$
$$0 = (\rho_m + k_m T) \cdot T_m$$

The steady states of the simplified system are located at:

- For \(T\) and \(T_m\):
  $$\{ (0, 0), \left(-\frac{\rho_m}{k_m}, \frac{\rho}{k_m} \right) \}$$

For this type of growth model, a trivial steady state (SS) exists. The steady state will change as the parameters change within this model, illustrating how the system's behavior responds to parameter variations.

# Future Work
For future work, is important that research will deepen into the in vivo analysis of HIV-1 dynamics to gain a better understanding of the virus's behavior within the host. Additionally, it will be valuable to investigate co-infection dynamics, examining how HIV-1 interacts with other diseases such as tuberculosis or hepatitis. Understanding these interactions is important for improving treatment protocols and patient outcomes. This broader scope of research will not only enhance the knowledge of HIV-1 but also inform the medical community on the best approaches to tackle the complexities of co-infections.

# Conclusions

In conclusion, the parameters estimated for the model align closely with those reported in the literature, providing a measure of validation. However, there remains a need to address the model's sensitivity to initial parameter guesses to ensure robustness. 

Sensitivity analysis has demonstrated that cell growth rate is the most critical parameter influencing the system's dynamics, underscoring its significance in the model. 

Bifurcation analysis revealed that the system exhibits two steady states. The occurrence of these states is contingent upon the accuracy of the model's parameter estimations, highlighting the importance of precise parameter determination. This understanding is instrumental in refining the model for future applications and studies. Future work should account for the analysis of the steady states without the complete initial set of ODEs of the system.

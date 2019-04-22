# Designing an Typhoid Environmental Surveillance and Optimize Sampling Site Allocation
Yuke Wang, Christine L. Moe, Shanta Dutta, Ashutosh Wadhwa, Suman Kanungo, Wolfgang Mairinger, Yichuan Zhao, Yi Jiang, and Peter F. M. Teunis

Environmental surveillance is used for monitoring disease in a population by detecting pathogens in sewage shedded by infected people. Detection of pathogens depends on many factors: infection rates and shedding in the population, and pathogen fate in the sewage network, but also sample size and sensitivity of assays. This complexity makes selection of sampling strategies challenging.

In the present study, a model was developed to simulate pathogen shedding, pathogen fate in the sewage network, environmental sampling, and detection of pathogens. The simulation study used <i>Salmonella enterica</i> serovar Typhi (<i>S.</i> Typhi) as the target pathogen and two wards in Kolkata, India as the study area.  Five different sampling strategies were evaluated for their sensitivities of detecting <i>S.</i> Typhi, by sampling units: pumping station, toilets, primary sampling units, pumping station + toilets, pumping station + primary sampling units. Sampling strategies were studied in 8 scenarios with different geographic clustering of risk, pathogen loss (decay, leakage), and sensitivity of detection assays. A novel adaptive sampling site allocation method was designed, which periodically updates the locations of sampling sites based on their performances. It is shown how the simulation model can be used to predict the performance of environmental surveillance and how it is improved by optimal allocation of sample sites. The simulation model has a number of settings that can be optimized to rapidly achieve high sensitivity: the initial number of sites, numbers of sites relocated per update, numbers of samples collected between updates of the proposed method can be optimized to rapidly achieve a sensitivity as high as possible.

The results were summarized as a decision tree to guide the selection of sampling strategy based on disease incidence, geographic distribution of risk, pathogen loss, and sensitivity of detection assays. The adaptive sampling site allocation method consistently outperformed alternatives with fixed site locations in most scenarios. In some cases, the optimum allocation method increased the median sensitivity from 45% to 90% within 20 updates.

<center><img src="cover.png" alt="Fecal" width="600"><center>

This repository includes all the codes (estimation, simulation) of Adaptive Sampling Site Allocation  model. The paper was submitted to Epidemics for publication.

## Acknowledgments
We are grateful to Anita Zaidi, Duncan Steele, Megan Carey, and Supriya Kumar for their insights in developing the SaniPath-Typhoid and Environmental Surveillance Strategy Study and the Bill & Melinda Gates Foundation for funding and supporting our work (grant OPP1150697). The study design, lab method development, and pilot data collection were conducted by the research team from the Center of Global Safe WASH, Emory University, including: Renuka Kapoor, Jamie Green, Suraja Raj, Pengbo Liu, Casey Siesel, Milagros Aldeco, Sarah Durry, and the research team from National Institute of Cholera and Enteric Diseases (NICED), including Ashish Mukhopadhyay, Pranab Chatterjee, Goutam Chowdhury. Sincere thanks to the Public Health Engineering Department, Govt. of West Bengal. for providing the sewage map for the city of Kolkata, India. We also want to thank John Crump for providing knowledge of Typhoid. We appreciate the help and insights from Lance Waller in building the mathematical simulation model of shedder dynamics.

## Questions
Please send any questions about working with Adaptive Sampling Site Allocation model to yuke.wang@emory.edu.

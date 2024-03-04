#ETBD
ETBD Simulaion 

## Setup
Requirements: R 3.0 or greater.




## Run simulations
The parameters governing a given simulation are specified here. These parameters include:

<table>
  <tr>
    <td>t</td>
    <td>number of time steps</td>
  </tr>
  <tr>
    <td>psymp</td>
    <td>probability of constant specation</td>
  </tr>
  <tr>
    <td>DIST</td>
    <td> Set a species abundnace distribution  log normal: "NORM" fishers log series: "SRS" geometric: "GEO none:"NO" </td>
  </tr>
  <tr>
    <td>watchgrow</td>
    <td>if TRUE allows animation of tree growth</td>
  </tr>
  <tr>
    <td>SADmarg</td>
    <td>the margin of error used in the SAD forcing mechanism</td>
  </tr>
  <tr>
    <td>siteN</td>
    <td>number of sites (needs to be the same number os the length of JmaxV) or a square if isGrid is TRUE)</td>
  </tr>
  <tr>
    <td>JmaxV</td>
    <td>vector of productivity zones (Jmax gradient) must be same length as numbe of sites</td>
  </tr>
   <tr>
    <td>NegExpEx</td>
    <td>if TRUE extinction is modeled as a negaitive exponential function</td>
  </tr>
    <tr>
    <td>exparm</td>
    <td>if NegExpEx is TRUE this is the negative exponential parameter for extinciton</td>
  </tr>
     <tr>
    <td>constantEX</td>
    <td>probability of extinction</td>
  </tr>
  <tr>
    <td>ExpSpParm</td>
    <td>if ExpSp is TRUE this is the exponential parameter for speciation</td>
  </tr>    
  <tr>
    <td>ExpSp</td>
    <td>if ExpSp is TRUE speciation is modeled as bounded exponential function</td>
  </tr>
  <tr>
   <td>splitparm</td>
    <td>the proportion of the population of the mother species that will split/or bud and be inheritated by the daughter species (controls the strength of heritability of population size)</td>
  </tr> 
  <tr>
   <td>split</td>
    <td>a small portion (depending on splitmarg) of the parent species will split off and migrate, default = TRUE</td>
  </tr>
  <tr>
     <td>bud</td>
    <td>if TRUE a small portion (depending on splitmarg) of the parent species will emerge leaving the origianl population of the parent species the same</td>
  </tr>
</table>








##simulation output 


The output files produced are listed below.

<table>
  <tr>
    <td>SIMULATION OUTPUT:</td>
    <td>a list of results</td>
  </tr>
  <tr>
    <td></td>
    <td>tree:   final tree</td>
  </tr>
  <tr>
    <td></td>
    <td>trees:   collection of every tree for time step</td>
  </tr>
  <tr>
    <td></td>
    <td>matrix_list:   final matrix of species locations and abundanecs</td>
  </tr>
  <tr>
    <td></td>
    <td>mig:   collection of matrix lists for every time step</td>  
  </tr>
  <tr>
    <td></td>
    <td>extincts:   list of all species that have gone extinct.</td>  
  </tr>
  <tr>
    <td></td>
    <td>allop:   matrix list of all species that allopatrically speciated</td>  
  </tr>
  <tr>
    <td></td>
    <td>symp:   matrix list of all species that have sympatrically speciated </td>  
  </tr>
  <tr>
    <td></td>
    <td>allotrip:   list of every allopatrically speciating species per time step</td>  
  </tr>
  <tr>
    <td></td>
    <td>symptrip:   list of every sympatrically speciating species per time step</td>  
  </tr>
  <tr>
    <td></td>
    <td>exsp:    matrix lists of extinct species per site </td>
  </tr>
</table>

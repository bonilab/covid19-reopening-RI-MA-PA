#### Format of output columns

The model has nine age classes, so output columns are groups by nines.  Group 1 is 0-9yo, Group 2 is 10-19yo, etc.

Some compartments will have multiple stages.  For example, the exposed class (E) has 6 stages in order to keep the variance in the incubation period low.  Thus, the exposed class is composed of 9 x 6 = 54 total population classes.



|   Columns    | Description | 
| -------------  |  ------------- |
| 1      | day (where day 1 is Jan 1 2020) |
| 2-10   | susceptible individuals |
| 11-64  | exposed individuals |
| 65-100 | asymptomatic individuals (also sub-clinical) |
| 101-136 | infected, symptomatic, infectious individuals -- but not hospitalized, or not yet hospitalized |
| 137-172 | hospitalized individuals, medical floor level of care (not ICU), in acute phase; i.e. individuals who were recently hospitalized; not individuals who have been discharged from the ICU because of improved condition |
| 173-181 | individuals in ICU but not on ventilator; acute phase; excludes individuals whose breathing has improved and who have been removed from mechanical ventilation but are still in the ICU |
| 182-235 | individuals in ICU and on ventilator ("invasive" mechanical ventilation or ECMO) |
| 236-244 | individuals in ICU who have recently been moved off a ventilator |
| 245-253 | hospitalized individuals, medical floor level of care (not ICU), in recovering phase; i.e. individuals whose condition in the ICU improved and have been moved outside of critical care |
| 254-262 | deceased individuals who died at home |
| 263-271 | deceased individuals who died in a hospital |
| 272-280 | recovered individuals who recovered from a non-hospitalized case of COVID19 (either symptomatic or asymptomatic) |
| 281-289 | recovered individuals who were discharged from a hospital and were hospitalized as a result of their COVID19 infection |
| 290-298 | cumulative symptomatic incidence by age |
| 299-307 | cumulative hospitalization incidence by age |

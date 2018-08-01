# PumaCoding2018
## Proyecto: obs_scheduling
#### Objetivo general:
Crear un script que comience varias observaciones consecutivas en forma automatica.
#### Estructura del script:
1. El observador ingresa el numero de observaciones N_obs .
2. El observador especifica, para cada observacion, los siguientes parametos del archivo de observacion \<pulsar\>.iar:
  * Highest Observation Frequency (MHz)
  * Observing time (minutes)
  * Average Data
3. Crear directorio del dia y subdirectorios para cada observacion (charlar convencion)
4. Comenzar la primer observacion. Al finalizar, moverla al subdirectorio correspondiente. Comenzar la segunda. Moverla. Etc.

**Extra:** Diseniar sistema de control para informar al observador de observaciones fallidas, atrasos, etc. Puede ser en pantalla, o a traves de un archivo *obs.info*.
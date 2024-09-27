# signalp-nextflow

Run signalp (https://services.healthtech.dtu.dk/services/SignalP-6.0/) versions 4,5,6 and combine the results.

Because version 6 is compute intensive we first filter the proteins by the signalp5 results.  The user can specify 
the minimum percentage of input proteins and the minimum SP score threshold.

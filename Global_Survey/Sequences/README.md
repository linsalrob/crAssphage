# Sequences

These are the raw DNA sequences provided by our [collaborators](../COLLABORATORS.md) worldwide.

The sequences are generated from three primer sequences that we designed, and using the general protocol below. Some people made minor variations to these protocols as they saw fit (e.g. replacing reagents).



## Primer sequences

* [Primer A](PrimerA):
  * Fwd: CTGATAGTATGATTGGTAAT
  * Rev: ATAAGTTCTCCAACTATCTT

* [Primer B](PrimerB):
  * Fwd: CCAGTATCTCCATAAGCATC
  * Rev: GTGAGGGCGGAATAGCTA

* [Primer C](PrimerC):
  * Fwd: GCAACAGGAGTAGTAAAATCTC
  * Rev: GCTCCTGTTAATCCTGATGTTA



To prepare the DNA template, we take raw sewage influent (the stuff coming into the sewage plant), centrifuge it briefly to remove the solids, and pass it through a 0.2 or 0.22 μm filter. We just use that in the reaction with no further treatment.

## PCR

The PCR mixture has

Reagent | Volume (μl)
--- | ---:
DNA template | 7.0
2x Master Mix |  | 25.0 μl
Forward primer (10μm) | 2.0 μl
Reverse primer (10μm) | 2.0 μl
DNAse free water  | 14.0 μl
Total volume | 50.0 μl


## Amplification

Amplification protocols:

Primer A:

* Denture 95°C for 3 minutes.

Then 30 cycles of:
  1. Denature 95°C for 45 seconds
  2. Annealing 42.6°C for 30 seconds
  3. Extension 68°C for 90 seconds

* Final extension	68°C for 5 minutes


Primers B & C:
* Denture 95°C for 3 minutes.

Then 30 cycles of:
  1. Denature 95°C for 45 seconds
  2. Annealing 50°C for 30 seconds
  3. Extension 68°C for 90 seconds

* Final extension	68°C for 5 minutes


**The expected product sizes are  A=1331 bp; B=1354 bp; C=1238 bp.**


# Metadata

We are trying to collect the metadata about where the sequences came from. Currently we are collecting the following information. Currently we include it in the header of each sequence. Although this results in a lot of redundancy, it also means that you have the information. There is no particular order to the keys, and not all sequences are guaranteed to have all the keys. If a key is missing it is assumed that the data is unknown.

* date: Ideally, the date that the sample was collected, but we will also use the date that the sequence was received if necessary.
* address: The location of the collection site. This is mainly used to verify lat/lon and we understand that not all locations want that shared!
* altitude: The height of the collection site above sea level
* latitude: The latitude of the collection site in [decimal degrees](https://en.wikipedia.org/wiki/Decimal_degrees)
* longitude: The longitude of the collection site in [decimal degrees](https://en.wikipedia.org/wiki/Decimal_degrees)
* name: The name of the sample or sequence. This can be an arbitrary string.



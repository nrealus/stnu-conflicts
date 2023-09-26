# STNU dynamic controllability checking with conflict extraction

A Python implementation of algorithms 1 & 2 from [1].

Dynamic controllability (DC) of STNUs is most often simply checked as being satisfied or not.
However, in the non-DC case, it can be interesting to retrieve some "conflicts" providing a reason as to why the STNU is not DC.
These take the form of a set (conjunction) of linear inequalities. The interest in them comes from the possibility to use them
as conflicts (disjunctions) in a conflict-driven or conflict-directed framework or solver. As a matter of fact, the code in this repository was
developed as a first step to integrate STNU DC-checking in a planner based on a hybrid CP/SAT solver. Moreover, as in [2], this
code may serve as a building block for PSTN risk-aware / chance-constrained DC-checking, which could be integrated in the planner as well.

Requires the [heapdict](https://pypi.org/project/HeapDict/) Python package for priority queues allowing priority updates.

TODO:
- More documentation relating to STNU theory (which I need a deeper understanding of).

References:

- [1]: [Bhargava, N., Vaquero, T. S., & Williams, B. C. (2017, August). Faster Conflict Generation for Dynamic Controllability. In IJCAI (pp. 4280-4286)](http://mers-papers.csail.mit.edu/Conference/2017/IJCAI17_Bhargava/RelaxIDC.pdf)

- [2]: [Wang, A.J. Risk-Bounded Dynamic Scheduling of Temporal Plans (2022) (PhD Thesis)](https://dspace.mit.edu/handle/1721.1/147542)

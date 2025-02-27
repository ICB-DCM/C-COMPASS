########
Glossary
########

This glossary provides definitions for terms and acronyms used in the context
of C-COMPASS.

.. glossary::
    :sorted:

    CC
    class contribution
        The relative fraction of a protein or other molecule type that is
        located in a specific compartment after filtering out false positives
        and re-normalizing the :term:`profile`.
        :math:`CC \in [0, 1]`.

    profile : B
        Here, a set of protein amounts across different fractions or
        compartments.

    CA
    class abundance
        Median of protein amounts of a given class.

    CC0
        The :term:`class contributions <class contribution>` before filtering#
        out false positives. I.e., the raw neural network output.
        :math:`CC0 \in [0, 1]`.

    DS
        Distance score.

    relocalization
    RL
        Relative change of protein localization for a specific compartment
        between two conditions. I.e., the difference between the
        :term:`class contributions <class contribution>` of a given protein
        in the two conditions.
        :math:`RL \in [-1, 1]`.

    relocalization score
    RLS
        Overall change of protein localization across all compartments
        (0 means no relocalization, 1 means 50% relocalization,
        2 means full relocalization) between two conditions,
        regardless of their origin and destination.
        :math:`RLS \in [0, 2]`.

    nCC
        Normalized class contribution (nCC = :term:`CC` * :term:`CA`).

    TPA
        Total protein amount.

    CPA
        Compartment protein amount (CPA = :term:`CC` * :term:`TPA`)

    nCPA
    normalized compartment protein amount
        Relative protein amount on a specific compartment,
        normalized by changes in compartment abundance and protein expression
        levels.
        This value can be used for comparison between conditions only.
        nCPA = :term:`nCC` * :term:`TPA`.

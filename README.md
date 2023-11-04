# Hole Filling Go
Go Implementation of the Mesh Hole Filling Algorithm by Peter Liepa.

This code is based on the work done by [russelmann](https://github.com/russelmann) in [this repo](https://github.com/russelmann/hole-filling-liepa).


# Using the Code

Using the code is straightforward. Specify the path to your mesh in the input to `readObj`. In the `FillHoleLiepa` function, you can choose the triangle weight calculation method to be either "area" or "angle", based on papers by: [Barequet & Sharir](https://dl.acm.org/doi/10.1016/0167-8396%2894%2900011-G) and [Peter Liepa](https://diglib.eg.org/handle/10.2312/SGP.SGP03.200-206), respectively.

**Please note:** Singular vertices are **not supported** by the `FindBoundaryLoops` function. That means no vertex should be adjacent to more than two boundary edges.

# Results

Left: original mesh - Middle: area-based method (Barequet paper) - Right: angle-based method (Liepa paper)

![Hole filling results](/img/allresults.png)


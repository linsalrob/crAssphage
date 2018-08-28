# The effect of different levels of alpha

We use the alpha channel (essentially the transparency) to show the number of times a line is drawn. For example, for the red lines, the alpha represents the number of times there are two sequences between two sites that are most similar to each other.

In the original version, we used an alpha of 0.2 that gives us this figure:

![World map with alpha=0.2](PrimerA_map.alpha0.2.png)

In response to a reviewer's question about why some sites didn't appear to have lines, we experimented with different levels of alpha. Here we provide the same image but with the alpha channel varying from 0.1 to 1.

Here is the image with alpha = 0.1

![World map with alpha=0.1](PrimerA_map.alpha0.1.png)

Here is the image with alpha = 1.0

![World map with alpha=1.0](PrimerA_map.alpha1.png)

We ended up using alpha = 0.3:

![World map with alpha=0.3](PrimerA_map.alpha0.3.png)




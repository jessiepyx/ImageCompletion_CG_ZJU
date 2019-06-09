### Pre-process of the dataset

##### 1. Mix the mask and origin image

(这里加两张图片：带mask的汽车，和origin汽车图片)

##### 2. Find the Adjacent Domain

- Traversing the image to generate a Map

  ![1560051175606](.\1560051175606.png)

### Region-filling algorithm

![img](./1560050536223.png?lastModify=1560050819?lastModify=1560050819)

##### 1. Find the edge points set

- Find the 8 points near one pixel in a square

- OR-Equal

  ![1560051984437](.\1560051984437.png)  

##### 2. Computing patch priority

- Filling order is crucial to non-parametric texture synthesis
- Best-first filling algorithm

这里要用两页PPT介绍（第一页介绍本文方法，第二页介绍后人方法，图片之后补充上去）

##### 2. Propagation texture information

- In specific area
- In neighborhood

![1560060665906](.\1560060665906.png)

##### 3. Updating confidence value

### Optimization for Texture Propagation

##### 1. Expended the Area to Fix

- Use a novel map to represent the expended area of the specific mask
- As for each edge point, the neighborhood of it should also be fixed

图片（待生成）

##### 2. Narrow the Area to Reference

- For each point, judge if it is the one with property texture
- When generating adjacent domain, set a safe-distance to keep away with structure line
- In order to avoid the points near the structure line with unclear texture

##### 3. Photometirc Correction






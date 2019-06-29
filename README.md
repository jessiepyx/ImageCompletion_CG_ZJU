# ImageCompletion_CG_ZJU

Implementation of 2005 SIGGRAPH paper [J. Sun et al. Image completion with structure propagation.](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/siggraph05_0265_final.pdf)


Code skeleton comes from [here](http://www.cad.zju.edu.cn/home/gfzhang/course/computational-photography/proj2-completion/completion.html).

**Collaborators:**

- [Jessie Peng](https://github.com/jessiepyx)
- [Xiu](https://github.com/Hap-Hugh)

## Structure Propagation

Contributed by Jessie Peng

- [x] User interface
  - choose input image
  - draw mask
  - draw structure lines/curve
  
- [x] Structure propagation
  - single line
  - multiple lines
  - curve
  
- [x] Photometric correction (copied from [https://github.com/ruanjiayi](https://github.com/ruanjiayi/Image-Completion-with-Structure-Propagation))

### Files

*cmake-build-debug/sp_result*: output after structure propagation

*cmake-build-debug/ts_result*: output after texture synthesis

*cmake-build-debug/mask_structure*: mask of structures

### Evironment

- Language: C++14
- IDE: CLion 2018.3
- OpenCV 4.0.1

## Texture Synthesis

Contributed by Xiu

Demo video in [weibo](https://m.weibo.cn/status/4380715229967145?wm=3333_2001&from=1095193010&sourcetype=qq&featurecode=newtitle) 
(Watching with your phone will be better)

### Evironment

- Language: C++14
- IDE: Microsoft Visual Studio 2017
- OpenCV 3.3.0 (suitable for OpenCV 4 and Opencv 2)

### Reference

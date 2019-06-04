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
  - draw structures
  
- [x] Structure propagation
  - single line
  - multiple lines
  - curve
  
- [ ] Photometric correction

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

# Moment method for thin loop wire

This Fortran code aims to solve a thin loop wire, of which the loop size is about 100 to 400 times the wavelength *lambda*.  
本程序采用 *Fortran* 语言编写，用以计算细导线圆环的散射场分布。目前完成进度可以计算出圆环上的电流分布。进一步的验证还没有进行。  
主程序位于 RCS 目录下，名为 *origin.f90*。 
使用到的工具主要是 *共轭梯度法（CGM）* 和 *快速傅里叶变换（FFT）*

## Conjugate Gradient Method
In this solution, the CGM must be available under the situation of complex numbers.  
在进行计算的时候，我们采用的共轭梯度算法需要对复数做一定的适应。迭代时间主要取决于矩阵和向量积（Matrix Vector Product）的复杂度。  
由于我们的求解目标是一个细圆环，其 **阻抗矩阵** 是一个具有 *Toeplitz* 特性的矩阵，因此该矩阵与任意向量的内积可以使用快速傅里叶变换( *Fast Fourier Transform* )来实现。

## Fast Fourier Transform
使用 库利－图基(Cooley-Tukey) 快速傅里叶变换算法  
将原来的矩阵做傅里叶变换，对待乘向量也进行傅里叶变换。变换后的两个向量对应位置乘积后再做反变换即可得到目标的矩阵向量积。  
[参考网页-维基百科](https://en.wikipedia.org/wiki/Cooley-Tukey_FFT_algorithm)  

## 时间及空间复杂度分析  

# SAR-SignalProcessing

# 合成孔径雷达算法与实现

## Catalogue

* Chapter_One

* Chapter_Two

    * Figure2.2&nbsp;&thinsp; 包含数据扭曲和旋转的傅里叶变换对

    * Figure2.3&nbsp;&thinsp; 矩形和sinc函数的傅里叶便函

    * Figure2.7&nbsp;&thinsp; 以两个不同采样率对300Hz的正弦波采样来说明混叠现象

    * Figure2.8&nbsp;&thinsp; 采样引起的频谱平移（实信号）

    * Figure2.9&nbsp;&thinsp; 采样引起的频谱平移（复信号）

    * Figure2.11 不同$\beta$值的Kaiser窗形状

    * Figure2.12 不同Kaiser窗的展宽和峰值旁瓣比

    * Figure2.14 使用sinc函数插值的图

    * FIgure2.15 Kaiser窗加权后的sinc函数，$\beta=2.5$

* Chapter_Three
 
    * Figure3.1&nbsp;&thinsp; 线性调频脉冲的相位和频率

    * Figure3.2&nbsp;&thinsp; 线性调频脉的复频谱

    * Figure3.3&nbsp;&thinsp; 不同TBP值的离散傅里叶变换频谱变化

    * Figure3.4&nbsp;&thinsp; 过采样率$\alpha_{os}$在频谱中引起的能量间隙

    * Figure3.5&nbsp;&thinsp; 匹配滤波器输出的3dB分辨率的测量

    * Figure3.6&nbsp;&thinsp; 基带线性调频信号的匹配滤波
    
    * Figure3.7&nbsp;&thinsp; 存在噪声时基带线性调频信号的匹配滤波
    
    * Figure3.8&nbsp;&thinsp; 非基带线性调频信号的匹配滤波
    
    * Figure3.9&nbsp;&thinsp; 匹配滤波后的信号频谱
    
    * Figure3.10 Kaiser窗在时域和频域中的实现形式
    
    * Figure3.11 方式2生成的频域匹配滤波器
    
    * Figure3.12 方式3生成的频域匹配滤波器
    
    * Figure3.13 通过压缩目标的位置来说明基带信号的弃置区和TA值

    * ~~Figure3.14 当$\beta=2.5$时的IRW、PSLR、ISLR与QPE之间的关系(结果并不吻合)~~

    * Figure3.15 旁瓣位置不同，而脉冲响应相似时的情况

* Chapter_Four

* Chapter_Five

    * Figure5.3&nbsp;&thinsp; 离散脉冲对方位信号的采样造成的方位混叠

    * Figure5.4&nbsp;&thinsp; 由方位chirp信号混叠造成的方位模糊

    * Figure5.5&nbsp;&thinsp; 斜视角为零和非零时的多普勒中心

    * Figure5.12 距离徙动的线性分量和二次分量

    * Figure5.13 目标轨迹在方位时域和方位频域的变化趋势

    * Figure5.16 零斜视角下单个点目标的时域性质
    
    * Figure5.17 零斜视角并且为正扫频下单个点目标的方位频谱
    
    * Figure5.18 非零斜视角下单个点目标的时域性质
    
    * Figure5.19 非零斜视角并且为正扫频下单个点目标的方位频谱

* Chapter_Six

    * Figure6.3&nbsp;&thinsp; 小斜视角情况下的多点雷达原始仿真信号

    * Figure6.4&nbsp;&thinsp; 距离压缩后的仿真结果

    * Figure6.5&nbsp;&thinsp; 方位向快速傅里叶变换后的仿真结果

    * Figure6.8&nbsp;&thinsp; 距离徙动校正不精确时，由调制引入的成对回波

    * Figure6.9&nbsp;&thinsp; 距离徙动校正后的仿真结果

    * Figure6.12 方位压缩后的仿真结果

    * Figure6.16 大斜视角情况下的多点雷达原始仿真信号

    * Figure6.17 方式3下二次距离压缩的精确实现

### Notes:

* .mlx实时脚本文件为练习时使用!

* arrow.m是一个非常牛逼的画箭头函数！
# Math for GPS Sim 


Install [Purple Pi](https://github.com/nschloe/purple-pi) to Render Math

[![purple-pi](https://img.shields.io/badge/Rendered%20with-Purple%20Pi-bd00ff?style=flat-square)](https://github.com/nschloe/purple-pi?activate)



## Introduction

(1)
$$ p_i = [(x-x_i)^2 + (y-y_i)^2 + (z - z_i)^2] + b = \rho_i +b \qquad$$

$$ p_i = c^* (t_{rec} - t_{sat} ) $$

## The spherical-plane algorithm

(2)
$$ (p_i- b)^2 = (x-x_i)^2 + (y-y_i)^2 + (z-z_i)^2, i =1,...,n \qquad$$

$$ (p_j - b)^2 - (p_j - b)^2 = [(x-x_j)^2 + (y-y_j)^2 + (z-z_j)^2] - [(x-x_{j+1})^2 + (y-y_{j+1})^2 + (z-z_{j+1})^2]
$$

(3)
$$ 0 = \alpha_j + \beta_j y + \gamma_j z + \delta_j b + \epsilon_j, j = 1, ..., n-1 \qquad$$

$$ \alpha_j = 2(x_{j+1} - x_j), \beta_j = 2(y_{j+1} - y_j), \\
\gamma_j = 2(z_{j+1} - z_j), \delta_j = 2(p_{1} - p_{j+1}), \\
\quad\\
\epsilon_j = (x_j^2 + y_j^2 + z_j^2) - (x_{j+1}^2 + y_{j+1}^2 + z_{j+1}^2) - (p_j^2 - p_{j+1}^2) = (r_j^2 - r_{j-1}^2) + (p_{j+1}^2 - p_j^2)
$$

### Solutions with five or more peudorange measurements

(4)
$$ (H^T H) \bm{s} = H^T \bm{u} \qquad$$

$$ H = \begin{bmatrix}\alpha_1 & \beta_1 & \gamma_1  & \delta_1 \\
       \alpha_2 & \beta_2 & \gamma_2  & \delta_2  \\
       \vdots  & \vdots  & \vdots & \vdots \\
       \alpha_{n-1} & \beta_{n-1} & \gamma_{n-1}  & \delta_{n-1}  
      \end{bmatrix},
      \bm{u} = \begin{bmatrix} \epsilon_1 \\ \epsilon_2 \\ \vdots \\ \epsilon_{n-1} \end{bmatrix}
$$

### Solutions with four peudorange measurements

$$ \begin{bmatrix}
      \alpha_1 & \beta_1 & \gamma_1  \\
       \alpha_2 & \beta_2 & \gamma_2   \\
       \alpha_3 & \beta_3 & \gamma_3
      \end{bmatrix}
      \begin{bmatrix} x \\ y \\ z \end{bmatrix} =
      \begin{bmatrix} \delta_1 \\ \delta_2 \\ \delta_3 \end{bmatrix} b -
      \begin{bmatrix} \epsilon_1 \\ \epsilon_2 \\ \epsilon_3 \end{bmatrix}
$$
or

(5)
$$
     \begin{bmatrix} x \\ y \\ z \end{bmatrix} =
      \begin{bmatrix} B_1 \\ B_2 \\ B_3 \end{bmatrix} b +
      \begin{bmatrix} C_1 \\ C_2 \\ C_3 \end{bmatrix} \qquad
$$

(6)
$0 = a_2x^2 + a_1x + a_0 \qquad$

where

$a_2 = B_1^2 + B_2^2 + B_3^2 - 1$

$a_1  = 2(p_i + B_1(C_1 - x_i) + B_2 (C_2 - y_i) + B_3 (C_3 - z_i))$

$a_0 = (C_1 - x_i)^2 + (C_2 -y_i)^2 + (C_3 - z_i)^2 - p^2_i$

(7)
$b = \frac{-a_1 + \sqrt{a_1^2 - 4a_2 a_0}}{2a_2} \qquad$

## The hyperbolic-plane algorithm

(8)
$\rho_i = \rho_j + (p_i - p_j)=  \rho_j + d_{ij} \qquad$

(9)
$\rho_i^2 = \rho_j^2 + 2 \rho_j d_ij + d_{ij}^2 \qquad$

(10)
$\rho_i^2 - \rho_j^2 = -2 * [x(x(x_i - x_j) + y(y_i - y_j) + z(z_i - z_j))] + r_i^2 - r_j^2 \qquad$

(11)
$\rho_j = \frac{1}{2 d_{ij}} \{ -2 [x (x_i - x_j) + y(y_i - y_j) + z(z_i - z_j)] + r_i^2 - r_j^2 - d^2_{ij} \} \qquad$

(12)
$\rho_j = \frac{1}{2 d_{kj}} \{ -2 [x (x_k - x_j) + y(y_k - y_j) + z(z_k - z_j)] + r_k^2 - r_j^2 - d_{kj}^2 \} \qquad$

$\rho_j^2 = (x-x_j)^2 + (y-y_j)^2 + (z-z_j)^2$

(13)
$P_2 x^2 + P_1 x + Q_2 y^2 + Q_1 y + R_2 z^2 + R_1 z + U_1 xy + U_2 xz + U_3 yz + W = 0 \qquad$
where

$$ P_2 = (x_1 - x_j)^2 - d_{ij}, P_1 = -e_{ij}(x_i  - x_j) + 2x_j d_{ij}^2, \\
Q_2 = (y_1 - y_j)^2 - d_{ij}, Q_1 = -e_{ij}(y_i  - y_j) + 2y_j d_{ij}^2, \\
R_2 = (z_1 - z_j)^2 - d_{ij}, R_1 = -e_{ij}(z_i  - z_j) + 2z_j d_{ij}^2, \\
U_1 = 2(x_i - x_j)(y_i - y_j), U_2 = 2(x_i - x_j) (z_i - z_j), \\
W = \frac{e_{ij}^2}{4} - d_{ij}^2 r_j^2  \\
e_ ij = r_i^2 - r_j^2 - d_{ij}^2
$$
(14)
$a_{ikj} x + b_{ikj} y + c_{ikj} z = f_{ikj} \qquad$

$$
a_{ikj} = \frac{x_i - x_j}{d_{ij}} - \frac{x_k- x_j}{d_{kj}}, \\
b_{ikj} = \frac{y_i - x_j}{d_{ij}} - \frac{x_k- x_j}{d_{kj}}, \\
c_{ikj} = \frac{z_i - z_j}{d_{ij}} - \frac{z_k- z_j}{d_{kj}}, \\
f_{ikj} = \frac{e_{ij}}{d_{ij}} - \frac{e_{kj}}{d_{kj}}
$$

### Solution with five or more peudorange measurements

(15)
$(H^T H) \bm{r} = H^T \bm{q} \qquad$

$$ H = \begin{bmatrix} a_1 & b_1 & c_1  \\
       a_2 & b_2 & c_2    \\
       \vdots  & \vdots  & \vdots \\
       a_{n-2} & b_{n-2} & c_{n-2}  
      \end{bmatrix},
      \bm{q} = \begin{bmatrix} a_{n-2} \\ b_{n-2} \\ c_{n-2} \end{bmatrix}
$$

(16)
$$b = \frac{1}{n} \sum_{i=1}^{n-2} (p_i - \rho_i) \qquad$$

### Solution with four peudorange measurements

$$  \begin{bmatrix} b_1 & c_1  \\
       b_2 & c_2     \\  
      \end{bmatrix} \begin{bmatrix} y \\ z  \end{bmatrix}
      = - \begin{bmatrix} a_1 \\ a_2 \end{bmatrix} x + \begin{bmatrix} f_1 \\ f_2 \end{bmatrix}
$$
or

(17)
$$  \begin{bmatrix} y \\ z  \end{bmatrix}
      = \begin{bmatrix} \alpha_1 \\ \alpha_2 \end{bmatrix} x + \begin{bmatrix} f
      \beta_1 \\ \beta_2 \end{bmatrix} \qquad $$

(18)
$0 = \gamma_2 x^2 + \gamma_1 x + \gamma_0 \qquad$

$$
\gamma_2 = P_2 + Q_2 \alpha_1^2 + R_2 \alpha_2 \alpha_2^2 + U_1 \alpha_2 + U_3 \alpha_1 \alpha_2 \\
\gamma_1 = P_1 + Q_1 \alpha_1  + R_1 \alpha_2 + U_1 \beta_1 + U_2 \beta_2 + 2(Q_2 \alpha_1 \beta_1 + R_2 \alpha_2 \beta_2) + U_3(\alpha_1 \beta_2 + \alpha_2 \beta_1)\\
\gamma_0 = W + Q_2 \beta_1^2 + Q_1 \beta_1 + R_2 \beta_2^2 + R_1 \beta_2 + U_3 \beta_1 \beta_2
$$
(19)
$x = \frac{-\gamma_1^2 - \sqrt{4 \gamma_2 \gamma_0}}{2 \gamma_2} \qquad$

## Sources

[How to write Latex in Markdown](http://flennerhag.com/2017-01-14-latex/)
[Alternative algorithms for the {GPS} static positioning solution by Lundberg, John B.](https://www.sciencedirect.com/science/article/pii/S0096300399002192)

## Scripts

<script type="text/javascript"
        src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS_CHTML"></script>

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
tex2jax: {
inlineMath: [['$','$'], ['\\(','\\)']],
processEscapes: true},
jax: ["input/TeX","input/MathML","input/AsciiMath","output/CommonHTML"],
extensions: ["tex2jax.js","mml2jax.js","asciimath2jax.js","MathMenu.js","MathZoom.js","AssistiveMML.js", "[Contrib]/a11y/accessibility-menu.js"],
TeX: {
extensions: ["AMSmath.js","AMSsymbols.js","noErrors.js","noUndefined.js"],
equationNumbers: {
autoNumber: "AMS"
}
}
});
</script>

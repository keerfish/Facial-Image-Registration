# Facial-Image-Registration
Facial image registration is to find distributed correspondences between two facial images accurately. Solving the facial image registration problem is a challenging task. The main reason is that the problem is ill-posed. The anisotropic elasticity of tissues under the human skin and the fact that deformations of the skin under facial expressions are based on the anatomic structures under the skin are two further difficulties in solving facial image
registration.

## Problem Statement
The image registration is supposed to occur between a facial image with a neutral expression and one with an angry expression as follows:
<p align="center">
  <img src="https://github.com/keerfish/Facial-Image-Registration/blob/main/imgs/input.jpg" align="center" width="650px"/>
</p>

The complex tissues under the skin make the skin mechanically anisotropic. However, the facial muscular system is responsible for the various facial expressions, especially the direction in which the muscle fibers dominate the movements from one facial expression to another.

## Main idea

Instead of using the strain energy of the anisotropic elastic material as a regularizing term, the strain energy of an isotropic elastic material can be projected in a pre-defined vector field. The pre-defined vector field characterizes the facial movement directions. The vector field may come from tracking the motion of some marked points on the face during a facial expression. In the end, the so-called projected strain energy and the strain energy of the anisotropic elastic material will have the same form under certain conditions.

A synthetic example is doing image registration between a black disc and a black oval. Since the deformation between the disc and oval is non-rigid. The main registration components utilized here are a B-Spline transformation, the SSD similarity metric, the bending energy penalty term, and the anatomy-adaptive regularizing term. The deformation from the black solid disc to the oval is illustrated below.
<p align="center">
  <img src="https://github.com/keerfish/Facial-Image-Registration/blob/main/imgs/synthetic.jpg" align="center" width="650px"/>
</p>

## Implementation
Within the framework of [elastix](http://elastix.isi.uu.nl;), we chose a reference image that comes
from a face with a neutral expression, and a template image, which comes from the same
face with an exaggerated expression. Under the anatomy-adaptive regularizing term, the
transformation tends to follow the proper directions.[elastix](http://elastix.isi.uu.nl;) is based on the [ITK](https://www.itk-engineering.de/gesundheit/robotik/?gad_source=1&gclid=CjwKCAiA1eO7BhATEiwAm0Ee-IeUb_wk3ex-3PqmOxi2apt0lpIKNHlmqJnWK8TsS6s4NpAudrZ-iBoCsnIQAvD_BwE), so the code base is thoroughly tested.  

When performing registration one carefully has to choose the several components. The
components must be specified in the parameter file. In the following table, a list of the components
that need to be specified is given, with some recommendations. The registration
component serves to connect all other components and implements the multiresolution aspect
of registration.

<p align="center">
  <img src="https://github.com/keerfish/Facial-Image-Registration/blob/main/imgs/components.jpg" align="center" width="600px"/>
</p>

## Results
The alignment of certain marked features serves as one of the criteria for evaluating the registration results. The green stars in the target image correspond to the landmarks in the reference image, while the blue stars represent the transformed landmarks after image registration.

<p align="center">
  <img src="https://github.com/keerfish/Facial-Image-Registration/blob/main/imgs/vectorfield.jpg" align="center" width="600px"/>
  <img src="https://github.com/keerfish/Facial-Image-Registration/blob/main/imgs/landmarkss.jpg" align="center" width="600px"/>
</p>






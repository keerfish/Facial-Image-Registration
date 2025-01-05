# Facial-Image-Registration
Facial image registration is to find distributed correspondences between two facial images accurately. Solving the facial image registration problem is a challenging task. The main reason is that the problem is ill-posed. The anisotropic elasticity of tissues under the human skin and the fact that deformations of the skin under facial expressions are based on the anatomic structures under the skin are two further difficulties in solving facial image
registration.

The complex tissues under the skin make the skin mechanically anisotropic. But the facial
muscular system is responsible for the various facial expressions, especially the direction in
which the muscle fibers run dominate the movements from one facial expression to another.

Instead of using the strain energy of the anisotropic elastic material as a regularizing
term, the strain energy of an isotropic elastic material can be projected in a pre-defined
vector field. The pre-defined vector field characterizes the facial movement directions. The
vector field may come from tracking the motion of some marked points on the face during a
facial expression. In the end, the so-called projected strain energy and the strain energy of
the anisotropic elastic material will have the same form under certain conditions.

Within the framework of elastix, we chose a reference image that comes
from a face with a neutral expression, and a template image, which comes from the same
face with an exaggerating expression. Under the anatomy-adaptive regularizing term, the
transformation tends to follow the proper directions.





/*=========================================================================
 *
 *  Copyright Zuse Institute Berlin and contributors
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkTransformDirectionalStrainEnergyPenaltyTerm_hxx
#define __itkTransformDirectionalStrainEnergyPenaltyTerm_hxx

#include "itkTransformDirectionalStrainEnergyPenaltyTerm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImage.h"

namespace itk {

/**
 * ****************** Constructor *******************************
 */

    template<class TFixedImage, class TScalarType>
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::TransformDirectionalStrainEnergyPenaltyTerm() {
        /** Pre-setting parameter. */
        this->m_PreSetParameter = NumericTraits<ScalarType>::Zero;

        /** Weights. */
        this->m_GeneralStrainPartWeight = NumericTraits<ScalarType>::One;
        this->m_PrincipleStrainPartWeight = NumericTraits<ScalarType>::One;

        this->m_horizontalDirectionWeight = NumericTraits<ScalarType>::One;
        this->m_verticalDirectionWeight = NumericTraits<ScalarType>::One;

        /** Values. */
        this->m_StrainEnergyPenaltyTermValue = NumericTraits<MeasureType>::Zero;

        this->m_GeneralStrainPartValue1 = NumericTraits<MeasureType>::Zero;
        this->m_GeneralStrainPartValue2 = NumericTraits<MeasureType>::Zero;

        this->m_PrincipleStrainPartValue1 = NumericTraits<MeasureType>::Zero;
        this->m_PrincipleStrainPartValue2 = NumericTraits<MeasureType>::Zero;
        this->m_PrincipleStrainPartValue3 = NumericTraits<MeasureType>::Zero;

        /** Gradient magnitudes. */
        this->m_GeneralStrainPartGradientMagnitude1 = NumericTraits<MeasureType>::Zero;
        this->m_GeneralStrainPartGradientMagnitude2 = NumericTraits<MeasureType>::Zero;

        this->m_PrincipleStrainPartGradientMagnitude1 = NumericTraits<MeasureType>::Zero;
        this->m_PrincipleStrainPartGradientMagnitude2 = NumericTraits<MeasureType>::Zero;
        this->m_PrincipleStrainPartGradientMagnitude3 = NumericTraits<MeasureType>::Zero;

        /** Usage. */
        this->m_UseGeneralStrainPart1 = true;
        this->m_UseGeneralStrainPart2 = true;

        this->m_UsePrincipleStrainPart1 = true;
        this->m_UsePrincipleStrainPart2 = true;
        this->m_UsePrincipleStrainPart3 = true;

        this->m_CalculateGeneralStrainPart1 = true;
        this->m_CalculateGeneralStrainPart2 = true;

        this->m_CalculatePrincipleStrainPart1 = true;
        this->m_CalculatePrincipleStrainPart2 = true;
        this->m_CalculatePrincipleStrainPart3 = true;

        /** Initialize dilation. */
        this->m_DilationRadiusMultiplier = NumericTraits<CoordinateRepresentationType>::One;
        this->m_DilateRigidityImages = true;
        this->m_DilateVectorImages = true;

        /** Initialize rigidity images and their usage. */
        this->m_UseFixedRigidityImage = true;
        this->m_UseMovingRigidityImage = true;
        this->m_FixedRigidityImage = 0;
        this->m_MovingRigidityImage = 0;

        this->m_UseHorizontalDirectionVectorImage = true;
        this->m_UseVerticalDirectionVectorImage = true;
        this->m_HorizontalDirectionVectorImage = 0;
        this->m_VerticalDirectionVectorImage = 0;

        this->m_RigidityCoefficientImage = RigidityImageType::New();
        this->m_RigidityCoefficientImageIsFilled = false;

        this->m_RigidityHorizontalVectorImage = VectorImageType1::New();
        this->m_RigidityHorizontalVectorImageIsFilled = false;

        this->m_RigidityVerticalVectorImage = VectorImageType2::New();
        this->m_RigidityVerticalVectorImageIsFilled = false;

        /** Initialize dilation filter for the rigidity images. */

        this->m_FixedRigidityImageDilation.resize(FixedImageDimension);
        this->m_MovingRigidityImageDilation.resize(MovingImageDimension);

        for (unsigned int i = 0; i < FixedImageDimension; i++) {
            this->m_FixedRigidityImageDilation[i] = 0;
            this->m_MovingRigidityImageDilation[i] = 0;
        }

        /** Initialize dilated rigidity images. */
        this->m_FixedRigidityImageDilated = 0;
        this->m_MovingRigidityImageDilated = 0;

        /** We don't use an image sampler for this advanced metric. */
        this->SetUseImageSampler(false);

        this->m_BSplineTransform = NULL;

        /** Initialize dilation filter for the horizontal direction images. */

        this->m_HorizontalDirectionVectorImageDilation.resize(FixedImageDimension);
        this->m_VerticalDirectionVectorImageDilation.resize(FixedImageDimension);

        for (unsigned int i = 0; i < FixedImageDimension; i++) {
            this->m_HorizontalDirectionVectorImageDilation[i] = 0;
            this->m_VerticalDirectionVectorImageDilation[i] = 0;
        }

        /** Initialize dilated rigidity images. */
        this->m_HorizontalDirectionVectorImageDilated = 0;
        this->m_VerticalDirectionVectorImageDilated = 0;

        /** We don't use an image sampler for this advanced metric. */
        this->SetUseImageSampler(false);

        this->m_BSplineTransform = NULL;


    } // end Constructor


/**
 * *********************** CheckUseAndCalculationBooleans *****************************
 */

    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::CheckUseAndCalculationBooleans(void) {
        if (this->m_UseGeneralStrainPart1) {
            this->m_CalculateGeneralStrainPart1 = true;
        }

        if (this->m_UseGeneralStrainPart2) {
            this->m_CalculateGeneralStrainPart2 = true;
        }

        if (this->m_UsePrincipleStrainPart1) {
            this->m_CalculatePrincipleStrainPart1 = true;
        }

        if (this->m_UsePrincipleStrainPart2) {
            this->m_CalculatePrincipleStrainPart2 = true;
        }

        if (this->m_UsePrincipleStrainPart3) {
            this->m_CalculatePrincipleStrainPart3 = true;
        }

    } // end CheckUseAndCalculationBooleans()


/**
 * *********************** Initialize *****************************
 */

    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::Initialize(void) throw(ExceptionObject) {
        /** Call the initialize of the superclass. */
        this->Superclass::Initialize();

        /** Check if this transform is a B-spline transform. */
        typename BSplineTransformType::Pointer localBSplineTransform = 0;
        bool transformIsBSpline = this->CheckForBSplineTransform(localBSplineTransform);
        if (transformIsBSpline) { this->SetBSplineTransform(localBSplineTransform); }

        /** Set the B-spline transform to m_RigidityPenaltyTermMetric. */
        if (!transformIsBSpline) {
            itkExceptionMacro( << "ERROR: this metric expects a B-spline transform." );
        }

        RigidityImageRegionType region;
        region.SetSize(localBSplineTransform->GetGridRegion().GetSize());
        region.SetIndex(localBSplineTransform->GetGridRegion().GetIndex());

        this->m_RigidityCoefficientImage->SetRegions(region);
        this->m_RigidityCoefficientImage->SetSpacing(
                localBSplineTransform->GetGridSpacing());
        this->m_RigidityCoefficientImage->SetOrigin(
                localBSplineTransform->GetGridOrigin());
        this->m_RigidityCoefficientImage->SetDirection(
                localBSplineTransform->GetGridDirection());
        this->m_RigidityCoefficientImage->Allocate();

        VectorImageRegionType1 region1;
        region1.SetSize(localBSplineTransform->GetGridRegion().GetSize());
        region1.SetIndex(localBSplineTransform->GetGridRegion().GetIndex());

        this->m_RigidityHorizontalVectorImage->SetRegions(region1);
        this->m_RigidityHorizontalVectorImage->SetSpacing(
                localBSplineTransform->GetGridSpacing());
        this->m_RigidityHorizontalVectorImage->SetOrigin(
                localBSplineTransform->GetGridOrigin());
        this->m_RigidityHorizontalVectorImage->SetDirection(
                localBSplineTransform->GetGridDirection());
        this->m_RigidityHorizontalVectorImage->Allocate();

        VectorImageRegionType2 region2;
        region2.SetSize(localBSplineTransform->GetGridRegion().GetSize());
        region2.SetIndex(localBSplineTransform->GetGridRegion().GetIndex());

        this->m_RigidityVerticalVectorImage->SetRegions(region2);
        this->m_RigidityVerticalVectorImage->SetSpacing(
                localBSplineTransform->GetGridSpacing());
        this->m_RigidityVerticalVectorImage->SetOrigin(
                localBSplineTransform->GetGridOrigin());
        this->m_RigidityVerticalVectorImage->SetDirection(
                localBSplineTransform->GetGridDirection());
        this->m_RigidityVerticalVectorImage->Allocate();

        if (!this->m_UseFixedRigidityImage && !this->m_UseMovingRigidityImage) {
            /** Fill the rigidity coefficient image with ones. */
            this->m_RigidityCoefficientImage->FillBuffer(1.0);
        } else {
            this->DilateRigidityImages();
        }

        /** Reset the filling bool. */
        this->m_RigidityCoefficientImageIsFilled = false;

        if (!this->m_UseHorizontalDirectionVectorImage) {
            this->m_RigidityHorizontalVectorImage->FillBuffer(128.0);
        } else {
            this->DilateVectorImages1();
        }

//  this->m_RigidityHorizontalVectorImage->MultiplyBuffer(1.0);

        this->m_RigidityHorizontalVectorImageIsFilled = false;


        if (!this->m_UseVerticalDirectionVectorImage) {
            this->m_RigidityVerticalVectorImage->FillBuffer(128.0);
        } else {
            this->DilateVectorImages2();
        }


        this->m_RigidityVerticalVectorImageIsFilled = false;


    } // end Initialize()



/**
 * **************** DilateRigidityImages *****************
 */

    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::DilateRigidityImages(void) {
        /** Dilate m_FixedRigidityImage and m_MovingRigidityImage. */
        if (this->m_DilateRigidityImages) {
            /** Some declarations. */
            SERadiusType radius;
            std::vector <StructuringElementType> structuringElement(FixedImageDimension);

            /** Setup the pipeline. */
            if (this->m_UseFixedRigidityImage) {
                /** Create the dilation filters for the fixedRigidityImage. */
                for (unsigned int i = 0; i < FixedImageDimension; i++) {
                    this->m_FixedRigidityImageDilation[i] = DilateFilterType::New();
                }
                this->m_FixedRigidityImageDilation[0]->SetInput(this->m_FixedRigidityImage);
            }
            if (this->m_UseMovingRigidityImage) {
                /** Create the dilation filter for the movingRigidityImage. */
                for (unsigned int i = 0; i < FixedImageDimension; i++) {
                    this->m_MovingRigidityImageDilation[i] = DilateFilterType::New();
                }
                this->m_MovingRigidityImageDilation[0]->SetInput(this->m_MovingRigidityImage);
            }

            /** Get the B-spline grid spacing. */
            GridSpacingType gridSpacing;
            if (this->m_BSplineTransform.IsNotNull()) {
                gridSpacing = this->m_BSplineTransform->GetGridSpacing();
            }

            /** Set stuff for the separate dilation. */
            for (unsigned int i = 0; i < FixedImageDimension; i++) {
                /** Create the structuring element. */
                radius.Fill(0);
                radius.SetElement(i,
                                  static_cast< unsigned long >(
                                          this->m_DilationRadiusMultiplier
                                          * gridSpacing[i] ));

                structuringElement[i].SetRadius(radius);
                structuringElement[i].CreateStructuringElement();

                /** Set the kernel into all dilation filters.
       * The SetKernel() is implemented using a itkSetMacro, so a
       * this->Modified() is automatically called, which is important,
       * since this changes every time Initialize() is called (every resolution).
       */

                if (this->m_UseFixedRigidityImage) {
                    this->m_FixedRigidityImageDilation[i]->SetKernel(structuringElement[i]);
                }
                if (this->m_UseMovingRigidityImage) {
                    this->m_MovingRigidityImageDilation[i]->SetKernel(structuringElement[i]);
                }

                /** Connect the pipelines. */
                if (i > 0) {
                    if (this->m_UseFixedRigidityImage) {
                        this->m_FixedRigidityImageDilation[i]->SetInput(
                                this->m_FixedRigidityImageDilation[i - 1]->GetOutput());
                    }
                    if (this->m_UseMovingRigidityImage) {
                        this->m_MovingRigidityImageDilation[i]->SetInput(
                                this->m_MovingRigidityImageDilation[i - 1]->GetOutput());
                    }
                }
            } // end for loop

            /** Do the dilation for m_FixedRigidityImage. */
            if (this->m_UseFixedRigidityImage) {
                try {
                    this->m_FixedRigidityImageDilation[FixedImageDimension - 1]->Update();
                }
                catch (itk::ExceptionObject &excp) {
                    /** Add information to the exception. */
                    excp.SetLocation("TransformDirectionalStrainEnergyPenaltyTerm - Initialize()");
                    std::string err_str = excp.GetDescription();
                    err_str += "\nError while dilating m_FixedRigidityImage.\n";
                    excp.SetDescription(err_str);
                    /** Pass the exception to an higher level. */
                    throw excp;
                }
            }

            /** Do the dilation for m_MovingRigidityImage. */
            if (this->m_UseMovingRigidityImage) {
                try {
                    this->m_MovingRigidityImageDilation[MovingImageDimension - 1]->Update();
                }
                catch (itk::ExceptionObject &excp) {
                    /** Add information to the exception. */
                    excp.SetLocation("TransformDirectionalStrainEnergyPenaltyTerm - Initialize()");
                    std::string err_str = excp.GetDescription();
                    err_str += "\nError while dilating m_MovingRigidityImage.\n";
                    excp.SetDescription(err_str);
                    /** Pass the exception to an higher level. */
                    throw excp;
                }
            }

            /** Put the output of the dilation into some dilated images. */
            if (this->m_UseFixedRigidityImage) {
                this->m_FixedRigidityImageDilated
                        = this->m_FixedRigidityImageDilation[FixedImageDimension - 1]->GetOutput();
            }
            if (this->m_UseMovingRigidityImage) {
                this->m_MovingRigidityImageDilated
                        = this->m_MovingRigidityImageDilation[MovingImageDimension - 1]->GetOutput();
            }
        } // end if rigidity images should be dilated
        else {
            /** Copy the pointers of the undilated images to the dilated ones
     * if no dilation is needed.
     */
            if (this->m_UseFixedRigidityImage) {
                this->m_FixedRigidityImageDilated = this->m_FixedRigidityImage;
            }
            if (this->m_UseMovingRigidityImage) {
                this->m_MovingRigidityImageDilated = this->m_MovingRigidityImage;
            }

        } // end else if

    } // end DilateRigidityImages()


/* **************** Begin DilateHorizontalVectorImages ****************/

    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::DilateVectorImages1(void) {
        /** Dilate m_HorizontalDirectionVectorImage and m_VerticalDirectionVectorImage. */
        if (this->m_DilateVectorImages) {
            /** Some declarations. */
            SERadiusType radius1;
            std::vector <StructuringElementType1> structuringElement1(FixedImageDimension);

            /** Setup the pipeline. */
            if (this->m_UseHorizontalDirectionVectorImage) {
                /** Create the dilation filters for the HorizontalDirectionVectorImage. */
                for (unsigned int i = 0; i < FixedImageDimension; i++) {
                    this->m_HorizontalDirectionVectorImageDilation[i] = DilateFilterType1::New();
                }
                this->m_HorizontalDirectionVectorImageDilation[0]->SetInput(this->m_HorizontalDirectionVectorImage);
            }


            /** Get the B-spline grid spacing. */
            GridSpacingType gridSpacing;
            if (this->m_BSplineTransform.IsNotNull()) {
                gridSpacing = this->m_BSplineTransform->GetGridSpacing();
            }

            /** Set stuff for the separate dilation. */
            for (unsigned int i = 0; i < FixedImageDimension; i++) {
                /** Create the structuring element. */
                radius1.Fill(0);
                radius1.SetElement(i,
                                   static_cast< unsigned long >(
                                           this->m_DilationRadiusMultiplier
                                           * gridSpacing[i] ));

                structuringElement1[i].SetRadius(radius1);
                structuringElement1[i].CreateStructuringElement();

                if (this->m_UseHorizontalDirectionVectorImage) {
                    this->m_HorizontalDirectionVectorImageDilation[i]->SetKernel(structuringElement1[i]);
                }

                /** Connect the pipelines. */
                if (i > 0) {
                    if (this->m_UseHorizontalDirectionVectorImage) {
                        this->m_HorizontalDirectionVectorImageDilation[i]->SetInput(
                                this->m_HorizontalDirectionVectorImageDilation[i - 1]->GetOutput());
                    }

                }
            } // end for loop

            /** Do the dilation for m_HorizontalDirectionVectorImage. */
            if (this->m_UseHorizontalDirectionVectorImage) {
                try {
                    this->m_HorizontalDirectionVectorImageDilation[FixedImageDimension - 1]->Update();
                }
                catch (itk::ExceptionObject &excp) {
                    /** Add information to the exception. */
                    excp.SetLocation("TransformDirectionalStrainEnergyPenaltyTerm - Initialize()");
                    std::string err_str = excp.GetDescription();
                    err_str += "\nError while dilating m_HorizontalDirectionVectorImage.\n";
                    excp.SetDescription(err_str);
                    /** Pass the exception to an higher level. */
                    throw excp;
                }
            }


            /** Put the output of the dilation into some dilated images. */
            if (this->m_UseHorizontalDirectionVectorImage) {
                this->m_HorizontalDirectionVectorImageDilated
                        = this->m_HorizontalDirectionVectorImageDilation[FixedImageDimension - 1]->GetOutput();
            }

        } // end if vector images should be dilated
        else {
            /** Copy the pointers of the undilated images to the dilated ones
     * if no dilation is needed.
     */
            if (this->m_UseHorizontalDirectionVectorImage) {
                this->m_HorizontalDirectionVectorImageDilated = this->m_HorizontalDirectionVectorImage;
            }

        } // end else if

    } // end DilateHorizontalVectorImages()


/* **************** Begin DilateVerticalVectorImages ****************/


    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::DilateVectorImages2(void) {
        /** Dilate m_HorizontalDirectionVectorImage . */
        if (this->m_DilateVectorImages) {
            /** Some declarations. */
            SERadiusType radius2;
            std::vector <StructuringElementType2> structuringElement2(FixedImageDimension);

            /** Setup the pipeline. */
            if (this->m_UseVerticalDirectionVectorImage) {
                /** Create the dilation filters for the HorizontalDirectionVectorImage. */
                for (unsigned int i = 0; i < FixedImageDimension; i++) {
                    this->m_VerticalDirectionVectorImageDilation[i] = DilateFilterType2::New();
                }
                this->m_VerticalDirectionVectorImageDilation[0]->SetInput(this->m_VerticalDirectionVectorImage);
            }

            /** Get the B-spline grid spacing. */
            GridSpacingType gridSpacing;
            if (this->m_BSplineTransform.IsNotNull()) {
                gridSpacing = this->m_BSplineTransform->GetGridSpacing();
            }

            /** Set stuff for the separate dilation. */
            for (unsigned int i = 0; i < FixedImageDimension; i++) {
                /** Create the structuring element. */
                radius2.Fill(0);
                radius2.SetElement(i,
                                   static_cast< unsigned long >(
                                           this->m_DilationRadiusMultiplier
                                           * gridSpacing[i] ));

                structuringElement2[i].SetRadius(radius2);
                structuringElement2[i].CreateStructuringElement();

                if (this->m_UseVerticalDirectionVectorImage) {
                    this->m_VerticalDirectionVectorImageDilation[i]->SetKernel(structuringElement2[i]);
                }

                /** Connect the pipelines. */
                if (i > 0) {

                    if (this->m_UseVerticalDirectionVectorImage) {
                        this->m_VerticalDirectionVectorImageDilation[i]->SetInput(
                                this->m_VerticalDirectionVectorImageDilation[i - 1]->GetOutput());
                    }
                }
            } // end for loop

            /** Do the dilation for m_VerticalDirectionVectorImage. */
            if (this->m_UseVerticalDirectionVectorImage) {
                try {
                    this->m_VerticalDirectionVectorImageDilation[FixedImageDimension - 1]->Update();
                }
                catch (itk::ExceptionObject &excp) {
                    /** Add information to the exception. */
                    excp.SetLocation("TransformDirectionalStrainEnergyPenaltyTerm - Initialize()");
                    std::string err_str = excp.GetDescription();
                    err_str += "\nError while dilating m_VerticalDirectionVectoryImage.\n";
                    excp.SetDescription(err_str);
                    /** Pass the exception to an higher level. */
                    throw excp;
                }
            }

            /** Put the output of the dilation into some dilated images. */

            if (this->m_UseVerticalDirectionVectorImage) {
                this->m_VerticalDirectionVectorImageDilated
                        = this->m_VerticalDirectionVectorImageDilation[FixedImageDimension - 1]->GetOutput();
            }
        } // end if rigidity images should be dilated
        else {
            /** Copy the pointers of the undilated images to the dilated ones
     * if no dilation is needed.
     */
            if (this->m_UseVerticalDirectionVectorImage) {
                this->m_VerticalDirectionVectorImageDilated = this->m_VerticalDirectionVectorImage;
            }

        } // end else if

    } // end DilateVectorImages()

/**
 * *************** FillRigidityCoefficientImage ****************
 */

    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::FillRigidityCoefficientImage(const ParametersType &parameters) const {
        /** Sanity check. */
        if (!this->m_UseFixedRigidityImage && !this->m_UseMovingRigidityImage) {
            return;
        }

        /** The rigidity image only changes when it depends on the moving image. */
        if (!this->m_UseMovingRigidityImage && this->m_RigidityCoefficientImageIsFilled) {
            return;
        }

        /** Make sure that the transform is up to date. */
        this->m_Transform->SetParameters(parameters);

        /** Create and reset an iterator over m_RigidityCoefficientImage. */
        RigidityImageIteratorType it(this->m_RigidityCoefficientImage,
                                     this->m_RigidityCoefficientImage->GetLargestPossibleRegion());
        it.GoToBegin();

        /** Fill m_RigidityCoefficientImage. */
        RigidityPixelType fixedValue, movingValue, in;
        RigidityImagePointType point;
        point.Fill(0.0f);
        RigidityImageIndexType index1, index2;
        index1.Fill(0);
        index2.Fill(0);
        fixedValue = NumericTraits<RigidityPixelType>::Zero;
        movingValue = NumericTraits<RigidityPixelType>::Zero;
        in = NumericTraits<RigidityPixelType>::Zero;
        bool isInFixedImage = false;
        bool isInMovingImage = false;

        while (!it.IsAtEnd()) {
            /** Get current pixel in world coordinates. */
            this->m_RigidityCoefficientImage
                    ->TransformIndexToPhysicalPoint(it.GetIndex(), point);

            /** Get the corresponding indices in the fixed and moving RigidityImage's.
     * NOTE: Floating point index results are truncated to integers.
     */

            if (this->m_UseFixedRigidityImage) {
                isInFixedImage = this->m_FixedRigidityImageDilated
                        ->TransformPhysicalPointToIndex(point, index1);

            }
            if (this->m_UseMovingRigidityImage) {
                isInMovingImage = this->m_MovingRigidityImageDilated
                        ->TransformPhysicalPointToIndex(
                                //this->m_Transform->TransformPoint( point ), index2 );
                                this->m_BSplineTransform->TransformPoint(point), index2);
            }

            /** Get the values at those positions. */
            if (this->m_UseFixedRigidityImage) {
                if (isInFixedImage) {
                    fixedValue = this->m_FixedRigidityImageDilated->GetPixel(index1);
                } else {
                    fixedValue = 0.0;
                }
            }

            if (this->m_UseMovingRigidityImage) {
                if (isInMovingImage) {
                    movingValue = this->m_MovingRigidityImageDilated->GetPixel(index2);

                } else {

                    movingValue = m_PreSetParameter;

                }

            }

            /** Determine the maximum. */
            if (this->m_UseFixedRigidityImage && this->m_UseMovingRigidityImage) {
                in = (fixedValue > movingValue ? fixedValue : movingValue);
            } else if (this->m_UseFixedRigidityImage && !this->m_UseMovingRigidityImage) {
                in = fixedValue;
            } else if (!this->m_UseFixedRigidityImage && this->m_UseMovingRigidityImage) {
                in = movingValue;
            }
            /** else{} is not happening here, because we assume that one of them is true.
     * In our case we checked that in the derived class: elxMattesMIWRR.
     */

            it.Set(in);

            /** Increase iterator. */
            ++it;
        } // end while loop over rigidity coefficient image

        /** Remember that the rigidity coefficient image is filled. */
        this->m_RigidityCoefficientImageIsFilled = true;


    } // end FillRigidityCoefficientImage()



/***************** FillRigidityHorizontalVectorImage ******************/

    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::FillRigidityHorizontalVectorImage(const ParametersType &parameters) const {
        /** Sanity check. */
        if (!this->m_UseHorizontalDirectionVectorImage) {
            return;
        }

        typename NearestNeighborInterpolatorType1::Pointer nearestNeighborInterpolatorType1 =
                NearestNeighborInterpolatorType1::New();
        typename ResampleFilterType::Pointer resizeFilter1 =
                ResampleFilterType::New();

        /** Make sure that the transform is up to date. */
        this->m_Transform->SetParameters(parameters);

        /** Create and reset an iterator over m_RigidityCoefficientImage. */
        VectorImageIteratorType1 it(this->m_RigidityHorizontalVectorImage,
                                    this->m_RigidityHorizontalVectorImage->GetLargestPossibleRegion());
        it.GoToBegin();

        /** Fill m_RigidityHorizontalVectorImage. */
        VectorPixelType1 fixedValue1, in1;
        VectorImagePointType1 point;
        point.Fill(0.0f);
        VectorImageIndexType1 index3;

        index3.Fill(0);
        fixedValue1 = NumericTraits<VectorPixelType1>::Zero;
        in1 = NumericTraits<VectorPixelType1>::Zero;

        bool isInHorizontalDirectionImage = false;

        while (!it.IsAtEnd()) {
            /** Get current pixel in world coordinates. */
            this->m_RigidityHorizontalVectorImage
                    ->TransformIndexToPhysicalPoint(it.GetIndex(), point);

            /** Get the corresponding indices in the fixed and moving RigidityImage's.
     * NOTE: Floating point index results are truncated to integers.
     */

            if (this->m_UseHorizontalDirectionVectorImage) {
                isInHorizontalDirectionImage = this->m_HorizontalDirectionVectorImageDilated
                        ->TransformPhysicalPointToIndex(point, index3);

            }

            /** Get the values at those positions. */
            if (this->m_UseHorizontalDirectionVectorImage) {
                if (isInHorizontalDirectionImage) {
                    fixedValue1 = this->m_HorizontalDirectionVectorImageDilated->GetPixel(index3);
                } else {
                    fixedValue1 = 128.0;
                    resizeFilter1->SetInterpolator(nearestNeighborInterpolatorType1);

                }
            }

            in1 = fixedValue1;

            it.Set(in1);

            /** Increase iterator. */
            ++it;
        } // end while loop.

        /** Remember that the Horizontal coefficient image is filled. */
        this->m_RigidityHorizontalVectorImageIsFilled = true;


    } // end FillRigidityHorizontalVectorImage()

/***************** FillRigidityVerticalVectorImage ******************/

    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::FillRigidityVerticalVectorImage(const ParametersType &parameters) const {
        /** Sanity check. */
        if (!this->m_UseVerticalDirectionVectorImage) {
            return;
        }

        typedef NearestNeighborInterpolateImageFunction<
                VectorImageType2, double> NearestNeighborInterpolatorType2;
        typedef typename NearestNeighborInterpolatorType2::Pointer nearestNeighborInterpolatorType2;

        /** Make sure that the transform is up to date. */
        this->m_Transform->SetParameters(parameters);

        /** Create and reset an iterator over m_RigidityCoefficientImage. */
        VectorImageIteratorType2 it(this->m_RigidityVerticalVectorImage,
                                    this->m_RigidityVerticalVectorImage->GetLargestPossibleRegion());
        it.GoToBegin();

        /** Fill m_RigidityHorizontalVectorImage. */
        VectorPixelType2 fixedValue2, in2;
        VectorImagePointType2 point;
        point.Fill(0.0f);
        VectorImageIndexType2 index5;
        index5.Fill(0);
        fixedValue2 = NumericTraits<VectorPixelType2>::Zero;
        in2 = NumericTraits<VectorPixelType2>::Zero;

        bool isInVerticalDirectionImage = false;
        while (!it.IsAtEnd()) {
            /** Get current pixel in world coordinates. */
            this->m_RigidityVerticalVectorImage
                    ->TransformIndexToPhysicalPoint(it.GetIndex(), point);

            /** Get the corresponding indices in the fixed and moving RigidityImage's.
     * NOTE: Floating point index results are truncated to integers.
     */

            if (this->m_UseVerticalDirectionVectorImage) {
                isInVerticalDirectionImage = this->m_VerticalDirectionVectorImageDilated
                        ->TransformPhysicalPointToIndex(point, index5);

            }

            /** Get the values at those positions. */
            if (this->m_UseVerticalDirectionVectorImage) {
                if (isInVerticalDirectionImage) {
                    fixedValue2 = this->m_VerticalDirectionVectorImageDilated->GetPixel(index5);
                } else {
                    fixedValue2 = 128.0;
                }
            }

            in2 = fixedValue2;

            it.Set(in2);

            /** Increase iterator. */
            ++it;
        } // end while loop over rigidity coefficient image

        /** Remember that the Vertical coefficient image is filled. */
        this->m_RigidityVerticalVectorImageIsFilled = true;


    } // end FillRigidityVerticalVectorImage()

/**
 * *********************** GetValue *****************************
 */
    template<class TFixedImage, class TScalarType>
    typename TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>::MeasureType
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::GetValue(const ParametersType &parameters) const {

        /** Fill the rigidity image based on the current transform parameters. */
        this->FillRigidityCoefficientImage(parameters);

        this->FillRigidityHorizontalVectorImage(parameters);

        this->FillRigidityVerticalVectorImage(parameters);

        /** Set output values to zero. */
        this->m_StrainEnergyPenaltyTermValue = NumericTraits<MeasureType>::Zero;

        this->m_GeneralStrainPartValue1 = NumericTraits<MeasureType>::Zero;
        this->m_GeneralStrainPartValue2 = NumericTraits<MeasureType>::Zero;

        this->m_PrincipleStrainPartValue1 = NumericTraits<MeasureType>::Zero;
        this->m_PrincipleStrainPartValue2 = NumericTraits<MeasureType>::Zero;
        this->m_PrincipleStrainPartValue3 = NumericTraits<MeasureType>::Zero;
        /** Set the parameters in the transform.
   * In this function, also the coefficient images are created.
   */
        this->m_BSplineTransform->SetParameters(parameters);

        /** Sanity check. */
        if (ImageDimension != 2) {
            itkExceptionMacro( << "ERROR: This filter is only implemented for dimension 2." );
        }

        /** Get a handle to the B-spline coefficient images. */
        std::vector <CoefficientImagePointer> inputImages(ImageDimension);
        for (unsigned int i = 0; i < ImageDimension; i++) {
            inputImages[i] = this->m_BSplineTransform->GetCoefficientImages()[i];
        }

        /** Get the B-spline coefficient image spacing. */
        CoefficientImageSpacingType spacing = inputImages[0]->GetSpacing();


        /** TASK 0:
   * Compute the rigidityCoefficientSum and check on it.
   *
   ************************************************************************* */

        /** Create iterator over the rigidity Coefficient image. */

        CoefficientImageIteratorType it_RCI(this->m_RigidityCoefficientImage,
                                            this->m_RigidityCoefficientImage->GetLargestPossibleRegion());
        it_RCI.GoToBegin();
        ScalarType rigidityCoefficientSum = NumericTraits<ScalarType>::Zero;


        HorizontalVectorImageIteratorType it_RCI1(this->m_RigidityHorizontalVectorImage,
                                                  this->m_RigidityHorizontalVectorImage->GetLargestPossibleRegion());
        it_RCI1.GoToBegin();
        ScalarType rigidityHorizontalVectorSum1 = NumericTraits<ScalarType>::Zero;


        VerticalVectorImageIteratorType it_RCI2(this->m_RigidityVerticalVectorImage,
                                                this->m_RigidityVerticalVectorImage->GetLargestPossibleRegion());
        it_RCI2.GoToBegin();
        ScalarType rigidityVerticalVectorSum2 = NumericTraits<ScalarType>::Zero;


        /** Add the rigidity coefficients together. */
        while (!it_RCI.IsAtEnd()) {

            rigidityCoefficientSum += it_RCI.Get();
            ++it_RCI;
        }

        /** Check for early termination*/
        if (rigidityCoefficientSum < 1e-14) {
            this->m_StrainEnergyPenaltyTermValue = NumericTraits<MeasureType>::Zero;
            return this->m_StrainEnergyPenaltyTermValue;
        }


        while (!it_RCI1.IsAtEnd()) {
            rigidityHorizontalVectorSum1 += it_RCI1.Get();
            ++it_RCI1;
        }

        /** Check for early termination. */
        if (rigidityHorizontalVectorSum1 < 1e-14) {
            this->m_StrainEnergyPenaltyTermValue = NumericTraits<MeasureType>::Zero;
            return this->m_StrainEnergyPenaltyTermValue;
        }


        while (!it_RCI2.IsAtEnd()) {
            rigidityVerticalVectorSum2 += it_RCI2.Get();
            ++it_RCI2;
        }

        /** Check for early termination. */
        if (rigidityVerticalVectorSum2 < 1e-14) {
            this->m_StrainEnergyPenaltyTermValue = NumericTraits<MeasureType>::Zero;
            return this->m_StrainEnergyPenaltyTermValue;
        }


        /** TASK 1:
   * Prepare for the calculation of the StrainEnergy penalty term.
   *
   ************************************************************************* */

        /** Create 1D neighborhood operators. */
        std::vector <NeighborhoodType> Operators_A(ImageDimension),
                Operators_B(ImageDimension),
                Operators_D(ImageDimension), Operators_E(ImageDimension),
                Operators_G(ImageDimension);

        /** Create B-spline coefficient images that are filtered once. */
        std::vector <CoefficientImagePointer> ui_FA(ImageDimension),
                ui_FB(ImageDimension), ui_FD(ImageDimension),
                ui_FE(ImageDimension), ui_FG(ImageDimension);

        /** For all dimensions ... */
        for (unsigned int i = 0; i < ImageDimension; i++) {
            /** ... create the filtered images ... */
            ui_FA[i] = CoefficientImageType::New();
            ui_FB[i] = CoefficientImageType::New();
            ui_FD[i] = CoefficientImageType::New();
            ui_FE[i] = CoefficientImageType::New();
            ui_FG[i] = CoefficientImageType::New();

            /** ... and the apropiate operators.
     * The operators C, D and E from the paper are here created
     * by Create1DOperator D, E and G, because of the 3D case and history.
     */
            this->Create1DOperator(Operators_A[i], "FA_xi", i + 1, spacing);
            this->Create1DOperator(Operators_B[i], "FB_xi", i + 1, spacing);
            this->Create1DOperator(Operators_D[i], "FD_xi", i + 1, spacing);
            this->Create1DOperator(Operators_E[i], "FE_xi", i + 1, spacing);
            this->Create1DOperator(Operators_G[i], "FG_xi", i + 1, spacing);

        } // end for loop

        /** TASK 2:
   * Filter the B-spline coefficient images.
   *
   ************************************************************************* */

        /** Filter the inputImages. */
        for (unsigned int i = 0; i < ImageDimension; i++) {
            ui_FA[i] = this->FilterSeparable(inputImages[i], Operators_A);
            ui_FB[i] = this->FilterSeparable(inputImages[i], Operators_B);
            ui_FD[i] = this->FilterSeparable(inputImages[i], Operators_D);
            ui_FE[i] = this->FilterSeparable(inputImages[i], Operators_E);
            ui_FG[i] = this->FilterSeparable(inputImages[i], Operators_G);

        }

        /** TASK 3:
   * Create iterators.
   *
   ************************************************************************* */

        /** Create iterators over ui_F?. */
        std::vector <CoefficientImageIteratorType> itA(ImageDimension),
                itB(ImageDimension), itD(ImageDimension),
                itE(ImageDimension), itG(ImageDimension);

        /** Create iterators. */
        for (unsigned int i = 0; i < ImageDimension; i++) {
            /** Create iterators. */
            itA[i] = CoefficientImageIteratorType(ui_FA[i], ui_FA[i]->GetLargestPossibleRegion());
            itB[i] = CoefficientImageIteratorType(ui_FB[i], ui_FB[i]->GetLargestPossibleRegion());
            itD[i] = CoefficientImageIteratorType(ui_FD[i], ui_FD[i]->GetLargestPossibleRegion());
            itE[i] = CoefficientImageIteratorType(ui_FE[i], ui_FE[i]->GetLargestPossibleRegion());
            itG[i] = CoefficientImageIteratorType(ui_FG[i], ui_FG[i]->GetLargestPossibleRegion());

            /** Reset iterators. */
            itA[i].GoToBegin();
            itB[i].GoToBegin();
            itD[i].GoToBegin();
            itE[i].GoToBegin();
            itG[i].GoToBegin();

        }

        /** TASK 4A:
   * Do the actual calculation of the StrainEnergy penalty term value.
   * Calculate the GeneralStrainPart term.
   *
   ************************************************************************* */

        /** Reset all iterators. */

        it_RCI.GoToBegin();
        it_RCI1.GoToBegin();
        it_RCI2.GoToBegin();

        if (this->m_CalculateGeneralStrainPart1) {
            ScalarType mu1_A, mu2_A, mu1_B, mu2_B;
            while (!itA[0].IsAtEnd()) {
                /** Copy values: this way we avoid calling Get() so many times.
       * It also improves code readability.
       */
                mu1_A = itA[0].Get();
                mu2_A = itA[1].Get();
                mu1_B = itB[0].Get();
                mu2_B = itB[1].Get();

                if (ImageDimension == 2) {
                    /** Calculate the value of the GeneralStrainPart1. */
                    this->m_GeneralStrainPartValue1
                            += it_RCI.Get() * (it_RCI1.Get() - 127) * (it_RCI1.Get() - 127)
                               * (
                                       vcl_pow(
                                               +(mu1_A / 2 + mu1_A / 2),
                                               2.0)
                                       + vcl_pow(
                                               +(mu2_A / 2 + mu1_B / 2),
                                               2.0)
                               );
                }


                /** Increase all iterators. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itA[i];
                    ++itB[i];

                }
                ++it_RCI;
                ++it_RCI1;
                ++it_RCI2;
            } //end while
        }// end if do GeneralStrainPart1

        it_RCI.GoToBegin();
        it_RCI2.GoToBegin();
        it_RCI1.GoToBegin();

        for (unsigned int i = 0; i < ImageDimension; i++) {
            itA[i].GoToBegin();
            itB[i].GoToBegin();

        }

        if (this->m_CalculateGeneralStrainPart2) {
            ScalarType mu1_A, mu2_A, mu1_B, mu2_B;
            while (!itA[0].IsAtEnd()) {
                /** Copy values: this way we avoid calling Get() so many times.
       * It also improves code readability.
       */
                mu1_A = itA[0].Get();
                mu2_A = itA[1].Get();
                mu1_B = itB[0].Get();
                mu2_B = itB[1].Get();


                if (ImageDimension == 2) {
                    /** Calculate the value of the GeneralStrainPart2. */
                    this->m_GeneralStrainPartValue2
                            += it_RCI.Get() * (it_RCI2.Get() - 127) * (it_RCI2.Get() - 127)
                               * (
                                       vcl_pow(
                                               +(mu1_B / 2 + mu2_A / 2),
                                               2.0)
                                       + vcl_pow(
                                               +(mu2_B + mu2_B),
                                               2.0)
                               );
                }


                /** Increase all iterators. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itA[i];
                    ++itB[i];

                }
                ++it_RCI;
                ++it_RCI2;
                ++it_RCI1;

            } // end while
        }   // end if do GeneralStrainPart2


        /** TASK 4B:
   * Do the actual calculation of the rigidity penalty term value.
   * Calculate the PrincipleStrain term.
   *
   ************************************************************************* */

        /** Reset all iterators. */
        for (unsigned int i = 0; i < ImageDimension; i++) {
            itA[i].GoToBegin();
            itB[i].GoToBegin();

        }
        it_RCI.GoToBegin();
        it_RCI1.GoToBegin();
        it_RCI2.GoToBegin();

        if (this->m_CalculatePrincipleStrainPart1) {
            ScalarType mu1_A, mu2_A, mu1_B, mu2_B;
            while (!itA[0].IsAtEnd()) {
                /** Copy values: this way we avoid calling Get() so many times.
       * It also improves code readability.
       */
                mu1_A = itA[0].Get();
                mu2_A = itA[1].Get();
                mu1_B = itB[0].Get();
                mu2_B = itB[1].Get();

                if (ImageDimension == 2) {
                    /** Calculate the value of the PrincipleStrainPart1. */
                    this->m_PrincipleStrainPartValue1
                            += it_RCI.Get() * (it_RCI1.Get() - 127) * (it_RCI1.Get() - 127) * (it_RCI1.Get() - 127) *
                               (it_RCI1.Get() - 127)
                               * (
                                       vcl_pow(
                                               +mu1_A, 2.0)
                               );
                }


                /** Increase all iterators. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itA[i];
                    ++itB[i];

                }
                ++it_RCI;
                ++it_RCI1;
                ++it_RCI2;

            } // end while
        }   // end if do PrincipleStrainPart1

        /** Reset all iterators. */
        for (unsigned int i = 0; i < ImageDimension; i++) {
            itA[i].GoToBegin();
            itB[i].GoToBegin();

        }
        it_RCI.GoToBegin();
        it_RCI2.GoToBegin();
        it_RCI1.GoToBegin();

        if (this->m_CalculatePrincipleStrainPart2) {
            ScalarType mu1_A, mu2_A, mu1_B, mu2_B;
            while (!itA[0].IsAtEnd()) {
                /** Copy values: this way we avoid calling Get() so many times.
         * It also improves code readability.
         */
                mu1_A = itA[0].Get();
                mu2_A = itA[1].Get();
                mu1_B = itB[0].Get();
                mu2_B = itB[1].Get();

                if (ImageDimension == 2) {

                    /** Calculate the value of the PrincipleStrainPart2. */
                    this->m_PrincipleStrainPartValue2
                            += it_RCI.Get() * (it_RCI2.Get() - 127) * (it_RCI2.Get() - 127) * (it_RCI2.Get() - 127) *
                               (it_RCI2.Get() - 127)
                               * (
                                       vcl_pow(
                                               +mu2_B, 2.0)
                               );
                }


                /** Increase all iterators. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itA[i];
                    ++itB[i];

                }
                ++it_RCI;
                ++it_RCI2;
                ++it_RCI1;

            } // end while
        }   // end if do PrincipleStrainPart2

        /** Reset all iterators. */
        for (unsigned int i = 0; i < ImageDimension; i++) {
            itA[i].GoToBegin();
            itB[i].GoToBegin();

        }
        it_RCI.GoToBegin();
        it_RCI1.GoToBegin();
        it_RCI2.GoToBegin();

        if (this->m_CalculatePrincipleStrainPart3) {
            ScalarType mu1_A, mu2_A, mu1_B, mu2_B, test, test1, test2;
            while (!itA[0].IsAtEnd()) {
                /** Copy values: this way we avoid calling Get() so many times.
           * It also improves code readability.
           */
                mu1_A = itA[0].Get();
                mu2_A = itA[1].Get();
                mu1_B = itB[0].Get();
                mu2_B = itB[1].Get();

                if (ImageDimension == 2) {
                    /** Calculate the value of the PrincipleStrainPart3. */

                    this->m_PrincipleStrainPartValue3
                            += it_RCI.Get() * (it_RCI1.Get() - 127) * (it_RCI2.Get() - 127) * (it_RCI1.Get() - 127) *
                               (it_RCI2.Get() - 127)
                               * (
                                       +(2 * mu1_A * mu2_B)
                               );

                }

                /** Increase all iterators. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itA[i];
                    ++itB[i];

                }
                ++it_RCI;
                ++it_RCI1;
                ++it_RCI2;

            } // end while
        }   // end if do PrincipleStrainPart3


        /** TASK 5:
   * Do the actual calculation of the StrainEnergy penalty term value.
   *
   ************************************************************************* */

        /** Calculate the StrainEnergy penalty term value. */
        if (this->m_CalculateGeneralStrainPart1) {
            this->m_GeneralStrainPartValue1 /= rigidityCoefficientSum;
        }

        if (this->m_CalculateGeneralStrainPart2) {
            this->m_GeneralStrainPartValue2 /= rigidityCoefficientSum;
        }

        if (this->m_CalculatePrincipleStrainPart1) {
            this->m_PrincipleStrainPartValue1 /= rigidityCoefficientSum;
        }

        if (this->m_CalculatePrincipleStrainPart2) {
            this->m_PrincipleStrainPartValue2 /= rigidityCoefficientSum;
        }

        if (this->m_CalculatePrincipleStrainPart3) {
            this->m_PrincipleStrainPartValue3 /= rigidityCoefficientSum;
        }


        if (this->m_UseGeneralStrainPart1) {
            this->m_StrainEnergyPenaltyTermValue
                    += this->m_GeneralStrainPartWeight * m_horizontalDirectionWeight
                       * this->m_GeneralStrainPartValue1;
        }

        if (this->m_UseGeneralStrainPart2) {
            this->m_StrainEnergyPenaltyTermValue
                    += this->m_GeneralStrainPartWeight * m_verticalDirectionWeight
                       * this->m_GeneralStrainPartValue2;
        }

        if (this->m_UsePrincipleStrainPart1) {
            this->m_StrainEnergyPenaltyTermValue
                    += this->m_PrincipleStrainPartWeight * m_horizontalDirectionWeight * m_horizontalDirectionWeight
                       * this->m_PrincipleStrainPartValue1;
        }

        if (this->m_UsePrincipleStrainPart2) {
            this->m_StrainEnergyPenaltyTermValue
                    += this->m_PrincipleStrainPartWeight * m_verticalDirectionWeight * m_verticalDirectionWeight
                       * this->m_PrincipleStrainPartValue2;
        }
        if (this->m_UsePrincipleStrainPart3) {
            this->m_StrainEnergyPenaltyTermValue
                    += this->m_PrincipleStrainPartWeight * m_horizontalDirectionWeight * m_verticalDirectionWeight
                       * this->m_PrincipleStrainPartValue3;
        }

//  std::cout << this->m_StrainEnergyPenaltyTermValue << std::endl;
        /** Return the rigidity penalty term value. */
        return this->m_StrainEnergyPenaltyTermValue;

    } // end GetValue()


/**
 * *********************** GetDerivative ************************
 */

    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::GetDerivative(const ParametersType &parameters,
                    DerivativeType &derivative) const {
        /** When the derivative is calculated, all information for calculating
   * the metric value is available. It does not cost anything to calculate
   * the metric value now. Therefore, we have chosen to only implement the
   * GetValueAndDerivative(), supplying it with a dummy value variable.
   */

        MeasureType dummyvalue = NumericTraits<MeasureType>::Zero;
        this->GetValueAndDerivative(parameters, dummyvalue, derivative);

    } // end GetDerivative()


/**
 * *********************** BeforeThreadedGetValueAndDerivative ***********************
 */

    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::BeforeThreadedGetValueAndDerivative(const TransformParametersType &parameters) const {
        /** In this function do all stuff that cannot be multi-threaded.
   * Meant for use in the combo-metric. So, I did not think about general usage yet.
   */
        if (this->m_UseMetricSingleThreaded) {
            this->m_BSplineTransform->SetParameters(parameters);
        }

    } // end BeforeThreadedGetValueAndDerivative()


/**
 * *********************** GetValueAndDerivative ****************
 */

    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::GetValueAndDerivative(const ParametersType &parameters,
                            MeasureType &value, DerivativeType &derivative) const {
        /** Fill the rigidity image based on the current transform parameters. */
        this->FillRigidityCoefficientImage(parameters);
        this->FillRigidityHorizontalVectorImage(parameters);
        this->FillRigidityVerticalVectorImage(parameters);
        /** Set output values to zero. */
        value = NumericTraits<MeasureType>::Zero;

        this->m_StrainEnergyPenaltyTermValue = NumericTraits<MeasureType>::Zero;

        this->m_GeneralStrainPartValue1 = NumericTraits<MeasureType>::Zero;
        this->m_GeneralStrainPartValue2 = NumericTraits<MeasureType>::Zero;

        this->m_PrincipleStrainPartValue1 = NumericTraits<MeasureType>::Zero;
        this->m_PrincipleStrainPartValue2 = NumericTraits<MeasureType>::Zero;
        this->m_PrincipleStrainPartValue3 = NumericTraits<MeasureType>::Zero;

        /** Set output values to zero. */
        derivative = DerivativeType(this->GetNumberOfParameters());
        derivative.Fill(NumericTraits<MeasureType>::ZeroValue());

        /** Call non-thread-safe stuff, such as:
   *   this->SetTransformParameters( parameters );
   *   this->GetImageSampler()->Update();
   * Because of these calls GetValueAndDerivative itself is not thread-safe,
   * so cannot be called multiple times simultaneously.
   * This is however needed in the CombinationImageToImageMetric.
   * In that case, you need to:
   * - switch the use of this function to on, using m_UseMetricSingleThreaded = true
   * - call BeforeThreadedGetValueAndDerivative once (single-threaded) before
   *   calling GetValueAndDerivative
   * - switch the use of this function to off, using m_UseMetricSingleThreaded = false
   * - Now you can call GetValueAndDerivative multi-threaded.
   */
        this->BeforeThreadedGetValueAndDerivative(parameters);

        /** Sanity check. */
        if (ImageDimension != 2) {
            itkExceptionMacro( << "ERROR: This filter is only implemented for dimension 2 ." );
        }

        /** Get a handle to the B-spline coefficient images. */
        std::vector <CoefficientImagePointer> inputImages(ImageDimension);
        for (unsigned int i = 0; i < ImageDimension; i++) {
            inputImages[i] = this->m_BSplineTransform->GetCoefficientImages()[i];
        }

        /** Get the B-spline coefficient image spacing. */
        CoefficientImageSpacingType spacing = inputImages[0]->GetSpacing();


        /** TASK 0:
   * Compute the rigidityCoefficientSum and check on it.
   *
   ************************************************************************* */

        /** Create iterator over the rigidity coefficient image. */


        CoefficientImageIteratorType it_RCI(this->m_RigidityCoefficientImage,
                                            this->m_RigidityCoefficientImage->GetLargestPossibleRegion());
        it_RCI.GoToBegin();
        ScalarType rigidityCoefficientSum = NumericTraits<ScalarType>::Zero;

        /** Add the rigidity coefficients together. */
        while (!it_RCI.IsAtEnd()) {
            rigidityCoefficientSum += it_RCI.Get();
            ++it_RCI;
        }

        /** Check for early termination. 不够严格，待优化 */
        if (rigidityCoefficientSum < 1e-14) {
            this->m_StrainEnergyPenaltyTermValue = NumericTraits<MeasureType>::Zero;
            return;
        }


        /** Create iterator over the rigidity coeficient image. */
        HorizontalVectorImageIteratorType it_RCI1(this->m_RigidityHorizontalVectorImage,
                                                  this->m_RigidityHorizontalVectorImage->GetLargestPossibleRegion());
        it_RCI1.GoToBegin();
        ScalarType rigidityHorizontalVectorSum1 = NumericTraits<ScalarType>::Zero;

        /** Add the rigidity coefficients together. */
        while (!it_RCI1.IsAtEnd()) {
            rigidityHorizontalVectorSum1 += it_RCI1.Get();
            ++it_RCI1;
        }

        /** Check for early termination. */
        if (rigidityHorizontalVectorSum1 < 1e-14) {
            this->m_StrainEnergyPenaltyTermValue = NumericTraits<MeasureType>::Zero;
            return;
        }


        /** Create iterator over the rigidity coefficient image. */
        VerticalVectorImageIteratorType it_RCI2(this->m_RigidityVerticalVectorImage,
                                                this->m_RigidityVerticalVectorImage->GetLargestPossibleRegion());
        it_RCI2.GoToBegin();
        ScalarType rigidityVerticalVectorSum2 = NumericTraits<ScalarType>::Zero;

        /** Add the rigidity coefficients together. */
        while (!it_RCI2.IsAtEnd()) {
            rigidityVerticalVectorSum2 += it_RCI2.Get();
            ++it_RCI2;
        }

        /** Check for early termination. */
        if (rigidityVerticalVectorSum2 < 1e-14) {
            this->m_StrainEnergyPenaltyTermValue = NumericTraits<MeasureType>::Zero;
            return;
        }



        /** TASK 1:
   * Prepare for the calculation of the rigidity penalty term.
   *
   ************************************************************************* */

        /** Create 1D neighborhood operators. */
        std::vector <NeighborhoodType> Operators_A(ImageDimension),
                Operators_B(ImageDimension), Operators_D(ImageDimension),
                Operators_E(ImageDimension), Operators_G(ImageDimension);

        /** Create B-spline coefficient images that are filtered once. */
        std::vector <CoefficientImagePointer> ui_FA(ImageDimension),
                ui_FB(ImageDimension), ui_FD(ImageDimension),
                ui_FE(ImageDimension), ui_FG(ImageDimension);

        /** For all dimensions ... */
        for (unsigned int i = 0; i < ImageDimension; i++) {
            /** ... create the filtered images ... */
            ui_FA[i] = CoefficientImageType::New();
            ui_FB[i] = CoefficientImageType::New();
            ui_FD[i] = CoefficientImageType::New();
            ui_FE[i] = CoefficientImageType::New();
            ui_FG[i] = CoefficientImageType::New();

            /** ... and the approiate operators.
     * The operators C, D and E from the paper are here created
     * by Create1DOperator D, E and G, because of the 3D case and history.
     */
            this->Create1DOperator(Operators_A[i], "FA_xi", i + 1, spacing);
            this->Create1DOperator(Operators_B[i], "FB_xi", i + 1, spacing);
            this->Create1DOperator(Operators_D[i], "FD_xi", i + 1, spacing);
            this->Create1DOperator(Operators_E[i], "FE_xi", i + 1, spacing);
            this->Create1DOperator(Operators_G[i], "FG_xi", i + 1, spacing);

        } // end for loop

        /** TASK 2:
   * Filter the B-spline coefficient images.
   *
   ************************************************************************* */

        /** Filter the inputImages. */
        for (unsigned int i = 0; i < ImageDimension; i++) {
            ui_FA[i] = this->FilterSeparable(inputImages[i], Operators_A);
            ui_FB[i] = this->FilterSeparable(inputImages[i], Operators_B);
            ui_FD[i] = this->FilterSeparable(inputImages[i], Operators_D);
            ui_FE[i] = this->FilterSeparable(inputImages[i], Operators_E);
            ui_FG[i] = this->FilterSeparable(inputImages[i], Operators_G);

        }

        /** TASK 3:
   * Create subparts and iterators.
   *
   ************************************************************************* */

        /** Create iterators over ui_F?. */
        std::vector <CoefficientImageIteratorType> itA(ImageDimension),
                itB(ImageDimension), itC(ImageDimension),
                itD(ImageDimension), itE(ImageDimension),
                itF(ImageDimension), itG(ImageDimension),
                itH(ImageDimension), itI(ImageDimension);

        /** Create iterators. */
        for (unsigned int i = 0; i < ImageDimension; i++) {
            /** Create iterators. */
            itA[i] = CoefficientImageIteratorType(ui_FA[i], ui_FA[i]->GetLargestPossibleRegion());
            itB[i] = CoefficientImageIteratorType(ui_FB[i], ui_FB[i]->GetLargestPossibleRegion());
            itD[i] = CoefficientImageIteratorType(ui_FD[i], ui_FD[i]->GetLargestPossibleRegion());
            itE[i] = CoefficientImageIteratorType(ui_FE[i], ui_FE[i]->GetLargestPossibleRegion());
            itG[i] = CoefficientImageIteratorType(ui_FG[i], ui_FG[i]->GetLargestPossibleRegion());

            /** Reset iterators. */
            itA[i].GoToBegin();
            itB[i].GoToBegin();
            itD[i].GoToBegin();
            itE[i].GoToBegin();
            itG[i].GoToBegin();

        }

        /** Create GeneralStrainPart and PrincipleStrainPart . */
        std::vector <std::vector<CoefficientImagePointer>> GSparts1(ImageDimension);
        std::vector <std::vector<CoefficientImagePointer>> GSparts2(ImageDimension);

        std::vector <std::vector<CoefficientImagePointer>> PSparts1(ImageDimension);
        std::vector <std::vector<CoefficientImagePointer>> PSparts2(ImageDimension);
        std::vector <std::vector<CoefficientImagePointer>> PSparts3(ImageDimension);

        for (unsigned int i = 0; i < ImageDimension; i++) {
            GSparts1[i].resize(ImageDimension);
            GSparts2[i].resize(ImageDimension);

            PSparts1[i].resize(ImageDimension);
            PSparts2[i].resize(ImageDimension);
            PSparts3[i].resize(ImageDimension);

            for (unsigned int j = 0; j < ImageDimension; j++) {
                GSparts1[i][j] = CoefficientImageType::New();
                GSparts1[i][j]->SetRegions(inputImages[0]->GetLargestPossibleRegion());
                GSparts1[i][j]->Allocate();

                GSparts2[i][j] = CoefficientImageType::New();
                GSparts2[i][j]->SetRegions(inputImages[0]->GetLargestPossibleRegion());
                GSparts2[i][j]->Allocate();

                PSparts1[i][j] = CoefficientImageType::New();
                PSparts1[i][j]->SetRegions(inputImages[0]->GetLargestPossibleRegion());
                PSparts1[i][j]->Allocate();

                PSparts2[i][j] = CoefficientImageType::New();
                PSparts2[i][j]->SetRegions(inputImages[0]->GetLargestPossibleRegion());
                PSparts2[i][j]->Allocate();

                PSparts3[i][j] = CoefficientImageType::New();
                PSparts3[i][j]->SetRegions(inputImages[0]->GetLargestPossibleRegion());
                PSparts3[i][j]->Allocate();

            }
        }


        /** Create iterators over all parts. */
        std::vector <std::vector<CoefficientImageIteratorType>> itGSp1(ImageDimension);
        std::vector <std::vector<CoefficientImageIteratorType>> itGSp2(ImageDimension);

        std::vector <std::vector<CoefficientImageIteratorType>> itPSp1(ImageDimension);
        std::vector <std::vector<CoefficientImageIteratorType>> itPSp2(ImageDimension);
        std::vector <std::vector<CoefficientImageIteratorType>> itPSp3(ImageDimension);

        for (unsigned int i = 0; i < ImageDimension; i++) {
            itGSp1[i].resize(ImageDimension);
            itGSp2[i].resize(ImageDimension);

            itPSp1[i].resize(ImageDimension);
            itPSp2[i].resize(ImageDimension);
            itPSp3[i].resize(ImageDimension);

            for (unsigned int j = 0; j < ImageDimension; j++) {
                itGSp1[i][j] = CoefficientImageIteratorType(GSparts1[i][j],
                                                            GSparts1[i][j]->GetLargestPossibleRegion());
                itGSp1[i][j].GoToBegin();

                itGSp2[i][j] = CoefficientImageIteratorType(GSparts2[i][j],
                                                            GSparts2[i][j]->GetLargestPossibleRegion());
                itGSp2[i][j].GoToBegin();

                itPSp1[i][j] = CoefficientImageIteratorType(PSparts1[i][j],
                                                            PSparts1[i][j]->GetLargestPossibleRegion());
                itPSp1[i][j].GoToBegin();

                itPSp2[i][j] = CoefficientImageIteratorType(PSparts2[i][j],
                                                            PSparts2[i][j]->GetLargestPossibleRegion());
                itPSp2[i][j].GoToBegin();

                itPSp3[i][j] = CoefficientImageIteratorType(PSparts3[i][j],
                                                            PSparts3[i][j]->GetLargestPossibleRegion());
                itPSp3[i][j].GoToBegin();
            }
        }

        /** TASK 4A:
   * Do the calculation of the GeneralStrain subparts.
   *
   ************************************************************************* */

        /** Reset all iterators. */
        it_RCI.GoToBegin();
        it_RCI1.GoToBegin();
        it_RCI2.GoToBegin();

        if (this->m_CalculateGeneralStrainPart1) {
            ScalarType mu1_A, mu2_A, mu1_B, mu2_B;
            ScalarType valueGS1;

            while (!itGSp1[0][0].IsAtEnd()) {
                /** Copy values: this way we avoid calling Get() so many times.
       * It also improves code readability.
       */
                mu1_A = itA[0].Get();
                mu2_A = itA[1].Get();
                mu1_B = itB[0].Get();
                mu2_B = itB[1].Get();

                if (ImageDimension == 2) {
                    /** Calculate the value of the GeneralStrainPart1. */
                    this->m_GeneralStrainPartValue1
                            += it_RCI.Get() * (it_RCI1.Get() - 127) * (it_RCI1.Get() - 127)
                               /*  /((it_RCI1.Get()-127)*(it_RCI1.Get()-127)+(it_RCI2.Get()-127)*(it_RCI2.Get()-127))*/
                               * (
                                       vcl_pow(
                                               +(mu1_A / 2 + mu1_A / 2),
                                               2.0)
                                       + vcl_pow(
                                               +(mu2_A / 2 + mu1_B / 2),
                                               2.0)
                               );

                    /** Calculate the derivative of the GeneralStrainPart. */

                    /** mu1, part A */

                    valueGS1
                            = +mu1_A;

                    itGSp1[0][0].Set(2.0 * valueGS1);

                    /** mu1, partB*/

                    valueGS1
                            = +(mu1_B / 2 + mu2_A / 2);

                    itGSp1[0][1].Set(valueGS1);

                    /** mu2, part A */
                    valueGS1
                            = +(mu1_B / 2 + mu2_A / 2);

                    itGSp1[1][0].Set(valueGS1);

                    /** mu2, partB*/
                    valueGS1
                            = +0;

                    itGSp1[1][1].Set(valueGS1);

                } // end if dim == 2

                /** Increase all iterators. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itA[i];
                    ++itB[i];

                    for (unsigned int j = 0; j < ImageDimension; j++) {
                        ++itGSp1[i][j];
                    }
                }
                ++it_RCI;
                ++it_RCI1;
                ++it_RCI2;

            } // end while
        }   // end if do GeneralStrainPart1

        /** Reset all iterators. */

        it_RCI.GoToBegin();
        it_RCI2.GoToBegin();
        it_RCI1.GoToBegin();

        for (unsigned int i = 0; i < ImageDimension; i++) {
            itA[i].GoToBegin();
            itB[i].GoToBegin();
        }


        if (this->m_CalculateGeneralStrainPart2) {
            ScalarType mu1_A, mu2_A, mu1_B, mu2_B;
            ScalarType valueGS2;

            while (!itGSp2[0][0].IsAtEnd()) {
                /** Copy values: this way we avoid calling Get() so many times.
         * It also improves code readability.
         */
                mu1_A = itA[0].Get();
                mu2_A = itA[1].Get();
                mu1_B = itB[0].Get();
                mu2_B = itB[1].Get();

                if (ImageDimension == 2) {
                    /** Calculate the value of the GeneralStrainPart2. */

                    this->m_GeneralStrainPartValue2
                            += it_RCI.Get() * (it_RCI2.Get() - 127) * (it_RCI2.Get() - 127)
                               /*  /((it_RCI1.Get()-127)*(it_RCI1.Get()-127)+(it_RCI2.Get()-127)*(it_RCI2.Get()-127)) */
                               * (
                                       vcl_pow(
                                               +(mu1_B / 2 + mu2_A / 2),
                                               2.0)
                                       + vcl_pow(
                                               +(mu2_B + mu2_B),
                                               2.0)
                               );

                    /** Calculate the derivative of the GeneralStrainPart. */

                    /** mu1, part A */

                    valueGS2
                            = +0;
                    itGSp2[0][0].Set(valueGS2);

                    /** mu1, partB*/

                    valueGS2
                            = +(mu1_B / 2 + mu2_A / 2);

                    itGSp2[0][1].Set(valueGS2);

                    /** mu2, part A */

                    valueGS2
                            = +(mu1_B / 2 + mu2_A / 2);

                    itGSp2[1][0].Set(valueGS2);

                    /** mu2, partB*/

                    valueGS2
                            = +mu2_B;

                    itGSp2[1][1].Set(valueGS2);

                } // end if dim == 2


                /** Increase all iterators. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itA[i];
                    ++itB[i];

                    for (unsigned int j = 0; j < ImageDimension; j++) {
                        ++itGSp2[i][j];
                    }
                }
                ++it_RCI;
                ++it_RCI2;
                ++it_RCI1;

            } // end while
        }   // end if do GeneralStrainPart2



        /** TASK 4B:
   * Do the calculation of the PrincipleStrainPart.
   *
   ************************************************************************* */

        /** Reset all iterators. */
        for (unsigned int i = 0; i < ImageDimension; i++) {
            itA[i].GoToBegin();
            itB[i].GoToBegin();
        }

        it_RCI.GoToBegin();
        it_RCI1.GoToBegin();
        it_RCI2.GoToBegin();

        if (this->m_CalculatePrincipleStrainPart1) {
            ScalarType mu1_A, mu2_A, mu1_B, mu2_B;
            ScalarType valuePS1;
            while (!itPSp1[0][0].IsAtEnd()) {
                /** Copy values: this way we avoid calling Get() so many times.
       * It also improves code readability.
       */
                mu1_A = itA[0].Get();
                mu2_A = itA[1].Get();
                mu1_B = itB[0].Get();
                mu2_B = itB[1].Get();

                if (ImageDimension == 2) {
                    /** Calculate the value of the PrincipleStrainPart. */
                    this->m_PrincipleStrainPartValue1
                            += it_RCI.Get() * (it_RCI1.Get() - 127) * (it_RCI1.Get() - 127) * (it_RCI1.Get() - 127) *
                               (it_RCI1.Get() - 127)
                               * (
                                       vcl_pow(
                                               +mu1_A, 2.0)
                               );

                    /** Calculate the derivative of the properness condition. 依次求偏导. */

                    /** mu1, part A*/
                    valuePS1
                            = +mu1_A;
                    itPSp1[0][0].Set(2.0 * valuePS1);

                    /** mu1, part B */
                    valuePS1
                            = +0;
                    itPSp1[0][1].Set(valuePS1);

                    /** mu2, part A */
                    valuePS1
                            = +0;
                    itPSp1[1][0].Set(valuePS1);
                    /** mu2, part B */
                    valuePS1
                            = +0;
                    itPSp1[1][1].Set(valuePS1);

                } // end if dim == 2


                /** Increase all iterators. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itA[i];
                    ++itB[i];

                    for (unsigned int j = 0; j < ImageDimension; j++) {
                        ++itPSp1[i][j];
                    }
                }
                ++it_RCI;
                ++it_RCI1;
                ++it_RCI2;

            } // end while
        }   // end if do PrincipleStrainPart1

        /** Reset all iterators. */
        for (unsigned int i = 0; i < ImageDimension; i++) {
            itA[i].GoToBegin();
            itB[i].GoToBegin();
        }

        it_RCI.GoToBegin();
        it_RCI2.GoToBegin();
        it_RCI1.GoToBegin();

        if (this->m_CalculatePrincipleStrainPart2) {
            ScalarType mu1_A, mu2_A, mu1_B, mu2_B, mu3_B;
            ScalarType valuePS2;
            while (!itPSp2[0][0].IsAtEnd()) {
                /** Copy values: this way we avoid calling Get() so many times.
         * It also improves code readability.
         */
                mu1_A = itA[0].Get();
                mu2_A = itA[1].Get();
                mu1_B = itB[0].Get();
                mu2_B = itB[1].Get();

                if (ImageDimension == 2) {
                    this->m_PrincipleStrainPartValue2
                            += it_RCI.Get() * (it_RCI2.Get() - 127) * (it_RCI2.Get() - 127) * (it_RCI2.Get() - 127) *
                               (it_RCI2.Get() - 127)
                               * (
                                       vcl_pow(
                                               +mu2_B, 2.0)
                               );

                    /** Calculate the derivative of the properness condition.*/

                    /** mu1, part A*/
                    valuePS2
                            = +0;
                    itPSp2[0][0].Set(valuePS2);

                    /** mu1, part B */
                    valuePS2
                            = +0;
                    itPSp2[0][1].Set(valuePS2);

                    /** mu2, part A */
                    valuePS2
                            = +0;
                    itPSp2[1][0].Set(valuePS2);
                    /** mu2, part B */
                    valuePS2
                            = +mu2_B;
                    itPSp2[1][1].Set(2.0 * valuePS2);

                } // end if dim == 2


                /** Increase all iterators. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itA[i];
                    ++itB[i];

                    for (unsigned int j = 0; j < ImageDimension; j++) {
                        ++itPSp2[i][j];
                    }
                }
                ++it_RCI;
                ++it_RCI2;
                ++it_RCI1;

            } // end while
        }   // end if do PrincipleStrainPart2

        /** Reset all iterators. */
        for (unsigned int i = 0; i < ImageDimension; i++) {
            itA[i].GoToBegin();
            itB[i].GoToBegin();
        }

        it_RCI.GoToBegin();
        it_RCI1.GoToBegin();
        it_RCI2.GoToBegin();

        if (this->m_CalculatePrincipleStrainPart3) {
            ScalarType mu1_A, mu2_A, mu1_B, mu2_B, tteesstt;
            ScalarType valuePS3;
            unsigned int i = 1;
            while (!itPSp3[0][0].IsAtEnd()) {
                /** Copy values: this way we avoid calling Get() so many times.
           * It also improves code readability.
           */
                mu1_A = itA[0].Get();
                mu2_A = itA[1].Get();
                mu1_B = itB[0].Get();
                mu2_B = itB[1].Get();

                if (ImageDimension == 2) {
                    /** Calculate the value of the PrincipleStrainPart. */

                    this->m_PrincipleStrainPartValue3
                            += it_RCI.Get() * (it_RCI1.Get() - 127) * (it_RCI2.Get() - 127) * (it_RCI1.Get() - 127) *
                               (it_RCI2.Get() - 127)
                               * (
                                       +(2 * mu1_A * mu2_B)
                               );

                    /** Calculate the derivative of the properness condition*/

                    /** mu1, part A*/
                    valuePS3
                            = +mu2_B;
                    itPSp3[0][0].Set(2.0 * valuePS3);

                    /** mu1, part B */
                    valuePS3
                            = +0;
                    itPSp3[0][1].Set(valuePS3);

                    /** mu2, part A */
                    valuePS3
                            = +0;
                    itPSp3[1][0].Set(valuePS3);
                    /** mu2, part B */
                    valuePS3
                            = +mu1_A;
                    itPSp3[1][1].Set(2.0 * valuePS3);

                } // end if dim == 2

                /** Increase all iterators. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itA[i];
                    ++itB[i];

                    for (unsigned int j = 0; j < ImageDimension; j++) {
                        ++itPSp3[i][j];
                    }
                }

                ++it_RCI;
                ++it_RCI1;
                ++it_RCI2;
                i++;

            } // end while
        }   // end if do PrincipleStrainPart3



        /** TASK 5:
   * Do the actual calculation of the rigidity penalty term value.
   *
   ************************************************************************* */

        /** Calculate the rigidity penalty term value. */

        if (this->m_CalculateGeneralStrainPart1) {
            this->m_GeneralStrainPartValue1 /= rigidityCoefficientSum;
        }

        if (this->m_CalculateGeneralStrainPart2) {
            this->m_GeneralStrainPartValue2 /= rigidityCoefficientSum;
        }

        if (this->m_CalculatePrincipleStrainPart1) {
            this->m_PrincipleStrainPartValue1 /= rigidityCoefficientSum;
        }

        if (this->m_CalculatePrincipleStrainPart2) {
            this->m_PrincipleStrainPartValue2 /= rigidityCoefficientSum;
        }

        if (this->m_CalculatePrincipleStrainPart3) {
            this->m_PrincipleStrainPartValue3 /= rigidityCoefficientSum;
        }

        if (this->m_UseGeneralStrainPart1) {
            this->m_StrainEnergyPenaltyTermValue
                    += this->m_GeneralStrainPartWeight * m_horizontalDirectionWeight
                       * this->m_GeneralStrainPartValue1;
        }

        if (this->m_UseGeneralStrainPart2) {
            this->m_StrainEnergyPenaltyTermValue
                    += this->m_GeneralStrainPartWeight * m_verticalDirectionWeight
                       * this->m_GeneralStrainPartValue2;
        }

        if (this->m_UsePrincipleStrainPart1) {
            this->m_StrainEnergyPenaltyTermValue
                    += this->m_PrincipleStrainPartWeight * m_horizontalDirectionWeight * m_horizontalDirectionWeight
                       * this->m_PrincipleStrainPartValue1;
        }

        if (this->m_UsePrincipleStrainPart2) {
            this->m_StrainEnergyPenaltyTermValue
                    += this->m_PrincipleStrainPartWeight * m_verticalDirectionWeight * m_verticalDirectionWeight
                       * this->m_PrincipleStrainPartValue2;
        }

        if (this->m_UsePrincipleStrainPart3) {
            this->m_StrainEnergyPenaltyTermValue
                    += this->m_PrincipleStrainPartWeight * m_horizontalDirectionWeight * m_verticalDirectionWeight
                       * this->m_PrincipleStrainPartValue3;
        }

        value = this->m_StrainEnergyPenaltyTermValue;

        /** TASK 6:
   * Create filtered versions of the subparts.
   * Create all necessary iterators and operators.
   ************************************************************************* */

        /** Create filtered GeneralStrainPart and PrincipleStrainPart. */
        std::vector <CoefficientImagePointer> GSpartsF1(ImageDimension);
        std::vector <CoefficientImagePointer> GSpartsF2(ImageDimension);

        std::vector <CoefficientImagePointer> PSpartsF1(ImageDimension);
        std::vector <CoefficientImagePointer> PSpartsF2(ImageDimension);
        std::vector <CoefficientImagePointer> PSpartsF3(ImageDimension);

        for (unsigned int i = 0; i < ImageDimension; i++) {
            GSpartsF1[i] = CoefficientImageType::New();
            GSpartsF1[i]->SetRegions(inputImages[0]->GetLargestPossibleRegion());
            GSpartsF1[i]->Allocate();

            GSpartsF2[i] = CoefficientImageType::New();
            GSpartsF2[i]->SetRegions(inputImages[0]->GetLargestPossibleRegion());
            GSpartsF2[i]->Allocate();

            PSpartsF1[i] = CoefficientImageType::New();
            PSpartsF1[i]->SetRegions(inputImages[0]->GetLargestPossibleRegion());
            PSpartsF1[i]->Allocate();

            PSpartsF2[i] = CoefficientImageType::New();
            PSpartsF2[i]->SetRegions(inputImages[0]->GetLargestPossibleRegion());
            PSpartsF2[i]->Allocate();

            PSpartsF3[i] = CoefficientImageType::New();
            PSpartsF3[i]->SetRegions(inputImages[0]->GetLargestPossibleRegion());
            PSpartsF3[i]->Allocate();
        }

        /** Create neighborhood iterators over the subparts. */
        std::vector <std::vector<NeighborhoodIteratorType>> nitGSp1(ImageDimension);
        std::vector <std::vector<NeighborhoodIteratorType>> nitGSp2(ImageDimension);

        std::vector <std::vector<NeighborhoodIteratorType>> nitPSp1(ImageDimension);
        std::vector <std::vector<NeighborhoodIteratorType>> nitPSp2(ImageDimension);
        std::vector <std::vector<NeighborhoodIteratorType>> nitPSp3(ImageDimension);

        RadiusType radius;
        radius.Fill(1);
        for (unsigned int i = 0; i < ImageDimension; i++) {
            nitGSp1[i].resize(ImageDimension);
            nitGSp2[i].resize(ImageDimension);

            nitPSp1[i].resize(ImageDimension);
            nitPSp2[i].resize(ImageDimension);
            nitPSp3[i].resize(ImageDimension);

            for (unsigned int j = 0; j < ImageDimension; j++) {
                nitGSp1[i][j] = NeighborhoodIteratorType(radius,
                                                         GSparts1[i][j], GSparts1[i][j]->GetLargestPossibleRegion());
                nitGSp1[i][j].GoToBegin();

                nitGSp2[i][j] = NeighborhoodIteratorType(radius,
                                                         GSparts2[i][j], GSparts2[i][j]->GetLargestPossibleRegion());
                nitGSp2[i][j].GoToBegin();

                nitPSp1[i][j] = NeighborhoodIteratorType(radius,
                                                         PSparts1[i][j], PSparts1[i][j]->GetLargestPossibleRegion());
                nitPSp1[i][j].GoToBegin();

                nitPSp2[i][j] = NeighborhoodIteratorType(radius,
                                                         PSparts2[i][j], PSparts2[i][j]->GetLargestPossibleRegion());
                nitPSp2[i][j].GoToBegin();

                nitPSp3[i][j] = NeighborhoodIteratorType(radius,
                                                         PSparts3[i][j], PSparts3[i][j]->GetLargestPossibleRegion());
                nitPSp3[i][j].GoToBegin();
            }
        }

        /** Create iterators over the filtered parts. */
        std::vector <CoefficientImageIteratorType> itGSpf1(ImageDimension);
        std::vector <CoefficientImageIteratorType> itGSpf2(ImageDimension);

        std::vector <CoefficientImageIteratorType> itPSpf1(ImageDimension);
        std::vector <CoefficientImageIteratorType> itPSpf2(ImageDimension);
        std::vector <CoefficientImageIteratorType> itPSpf3(ImageDimension);

        for (unsigned int i = 0; i < ImageDimension; i++) {
            itGSpf1[i] = CoefficientImageIteratorType(GSpartsF1[i],
                                                      GSpartsF1[i]->GetLargestPossibleRegion());
            itGSpf1[i].GoToBegin();

            itGSpf2[i] = CoefficientImageIteratorType(GSpartsF2[i],
                                                      GSpartsF2[i]->GetLargestPossibleRegion());
            itGSpf2[i].GoToBegin();

            itPSpf1[i] = CoefficientImageIteratorType(PSpartsF1[i],
                                                      PSpartsF1[i]->GetLargestPossibleRegion());
            itPSpf1[i].GoToBegin();

            itPSpf2[i] = CoefficientImageIteratorType(PSpartsF2[i],
                                                      PSpartsF2[i]->GetLargestPossibleRegion());
            itPSpf2[i].GoToBegin();

            itPSpf3[i] = CoefficientImageIteratorType(PSpartsF3[i],
                                                      PSpartsF3[i]->GetLargestPossibleRegion());
            itPSpf3[i].GoToBegin();

        }

        /** Create a neighborhood iterator over the rigidity image. */
        NeighborhoodIteratorType nit_RCI(radius, this->m_RigidityCoefficientImage,
                                         this->m_RigidityCoefficientImage->GetLargestPossibleRegion());
        nit_RCI.GoToBegin();
        unsigned int neighborhoodSize = nit_RCI.Size();

        NeighborhoodIteratorType1 nit_RCI1(radius, this->m_RigidityHorizontalVectorImage,
                                           this->m_RigidityHorizontalVectorImage->GetLargestPossibleRegion());
        nit_RCI1.GoToBegin();
        unsigned int neighborhoodSize1 = nit_RCI1.Size();

        NeighborhoodIteratorType2 nit_RCI2(radius, this->m_RigidityVerticalVectorImage,
                                           this->m_RigidityVerticalVectorImage->GetLargestPossibleRegion());
        nit_RCI2.GoToBegin();
        unsigned int neighborhoodSize2 = nit_RCI2.Size();

        /** Create ND operators. */
        NeighborhoodType Operator_A, Operator_B, Operator_D, Operator_E,
                Operator_G;
        this->CreateNDOperator(Operator_A, "FA", spacing);
        this->CreateNDOperator(Operator_B, "FB", spacing);



        /** TASK 7A:
   * Calculate the filtered versions of the GeneralStrainPart subparts.
   * These are F_A * {subpart_0} + F_B * {subpart_1},
   * and (for 3D) + F_C * {subpart_2}, for all dimensions.
   ************************************************************************* */

        if (this->m_CalculateGeneralStrainPart1) {
            while (!itGSpf1[0].IsAtEnd()) {
                /** Create and reset tmp with zeros. */
                std::vector<double> tmp(ImageDimension, 0.0);

                /** Loop over all dimensions. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    /** Loop over the neighborhood. */
                    for (unsigned int k = 0; k < neighborhoodSize; ++k) {
                        /** Calculation of the inner product. */
                        tmp[i] += Operator_A.GetElement(k)      // FA *
                                  * nitGSp1[i][0].GetPixel(k)          // subpart[ i ][ 0 ]
                                  * nit_RCI.GetPixel(k)
                                  * (nit_RCI1.GetPixel(k) - 127)
                                  * (nit_RCI1.GetPixel(k) - 127);
                        tmp[i] += Operator_B.GetElement(k)      // FB *
                                  * nitGSp1[i][1].GetPixel(k)          // subpart[ i ][ 1 ]
                                  * nit_RCI.GetPixel(k)
                                  * (nit_RCI1.GetPixel(k) - 127)
                                  * (nit_RCI1.GetPixel(k) - 127);

                    } // end loop over neighborhood

                    /** Set the result in the filtered part. */
                    itGSpf1[i].Set(tmp[i]);

                } // end loop over dimension i

                /** Increase all iterators. */
                ++nit_RCI;
                ++nit_RCI1;
                ++nit_RCI2;

                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itGSpf1[i];
                    for (unsigned int j = 0; j < ImageDimension; j++) {
                        ++nitGSp1[i][j];
                    }
                }
            } // end while
        }   // end if do GeneralStrainPart1

        nit_RCI.GoToBegin();
        nit_RCI1.GoToBegin();
        nit_RCI2.GoToBegin();

        if (this->m_CalculateGeneralStrainPart2) {
            while (!itGSpf2[0].IsAtEnd()) {
                /** Create and reset tmp with zeros. */
                std::vector<double> tmp(ImageDimension, 0.0);

                /** Loop over all dimensions. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    /** Loop over the neighborhood. */
                    for (unsigned int k = 0; k < neighborhoodSize; ++k) {
                        /** Calculation of the inner product. */
                        tmp[i] += Operator_A.GetElement(k)      // FA *
                                  * nitGSp2[i][0].GetPixel(k)          // subpart[ i ][ 0 ]
                                  * nit_RCI.GetPixel(k)
                                  * (nit_RCI2.GetPixel(k) - 127)
                                  * (nit_RCI2.GetPixel(k) - 127)
                            /*	/((nit_RCI1.GetPixel( k )-127)*(nit_RCI1.GetPixel( k )-127)+(nit_RCI2.GetPixel( k )-127)*(nit_RCI2.GetPixel( k )-127)) */;
                        // c(k)
                        tmp[i] += Operator_B.GetElement(k)      // FB *
                                  * nitGSp2[i][1].GetPixel(k)          // subpart[ i ][ 1 ]
                                  * nit_RCI.GetPixel(k)
                                  * (nit_RCI2.GetPixel(k) - 127)
                                  * (nit_RCI2.GetPixel(k) - 127);
                    } // end loop over neighborhood

                    /** Set the result in the filtered part. */
                    itGSpf2[i].Set(tmp[i]);

                } // end loop over dimension i

                /** Increase all iterators. */
                ++nit_RCI;
                ++nit_RCI2;
                ++nit_RCI1;
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itGSpf2[i];
                    for (unsigned int j = 0; j < ImageDimension; j++) {
                        ++nitGSp2[i][j];
                    }
                }
            } // end while
        }   // end if do GeneralStrainPart2

        /** TASK 7B:
   * Calculate the filtered versions of the PrincipleStrainPart subparts.
   * These are F_A * {subpart_0} + F_B * {subpart_1},
   * and (for 3D) + F_C * {subpart_2}, for all dimensions.
   ************************************************************************* */

        nit_RCI.GoToBegin();
        nit_RCI1.GoToBegin();
        nit_RCI2.GoToBegin();

        if (this->m_CalculatePrincipleStrainPart1) {
            while (!itPSpf1[0].IsAtEnd()) {
                /** Create and reset tmp with zeros. */
                std::vector<double> tmp(ImageDimension, 0.0);

                /** Loop over all dimensions. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    /** Loop over the neighborhood. */
                    for (unsigned int k = 0; k < neighborhoodSize; ++k) {
                        /** Calculation of the inner product. */
                        tmp[i] += Operator_A.GetElement(k)      // FA *
                                  * nitPSp1[i][0].GetPixel(k)          // subpart[ i ][ 0 ]
                                  * nit_RCI.GetPixel(k)
                                  * (nit_RCI1.GetPixel(k) - 127)
                                  * (nit_RCI1.GetPixel(k) - 127)
                                  * (nit_RCI1.GetPixel(k) - 127)
                                  * (nit_RCI1.GetPixel(k) - 127);

                        tmp[i] += Operator_B.GetElement(k)      // FB *
                                  * nitPSp1[i][1].GetPixel(k)          // subpart[ i ][ 1 ]
                                  * nit_RCI.GetPixel(k)
                                  * (nit_RCI1.GetPixel(k) - 127)
                                  * (nit_RCI1.GetPixel(k) - 127)
                                  * (nit_RCI1.GetPixel(k) - 127)
                                  * (nit_RCI1.GetPixel(k) - 127);

                    } // end loop over neighborhood

                    /** Set the result in the filtered part. */
                    itPSpf1[i].Set(tmp[i]);

                } // end loop over dimension i

                /** Increase all iterators. */
                ++nit_RCI;
                ++nit_RCI1;
                ++nit_RCI2;

                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itPSpf1[i];
                    for (unsigned int j = 0; j < ImageDimension; j++) {
                        ++nitPSp1[i][j];
                    }
                }
            } // end while
        }   // end if do PrincipleStrainPart1

        nit_RCI.GoToBegin();
        nit_RCI2.GoToBegin();
        nit_RCI1.GoToBegin();

        if (this->m_CalculatePrincipleStrainPart2) {
            while (!itPSpf2[0].IsAtEnd()) {
                /** Create and reset tmp with zeros. */
                std::vector<double> tmp(ImageDimension, 0.0);

                /** Loop over all dimensions. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    /** Loop over the neighborhood. */
                    for (unsigned int k = 0; k < neighborhoodSize; ++k) {
                        /** Calculation of the inner product. */
                        tmp[i] += Operator_A.GetElement(k)      // FA *
                                  * nitPSp2[i][0].GetPixel(k)          // subpart[ i ][ 0 ]
                                  * nit_RCI.GetPixel(k)
                                  * (nit_RCI2.GetPixel(k) - 127)
                                  * (nit_RCI2.GetPixel(k) - 127)
                                  * (nit_RCI2.GetPixel(k) - 127)
                                  * (nit_RCI2.GetPixel(k) - 127);

                        tmp[i] += Operator_B.GetElement(k)      // FB *
                                  * nitPSp2[i][1].GetPixel(k)          // subpart[ i ][ 1 ]
                                  * nit_RCI.GetPixel(k)
                                  * (nit_RCI2.GetPixel(k) - 127)
                                  * (nit_RCI2.GetPixel(k) - 127)
                                  * (nit_RCI2.GetPixel(k) - 127)
                                  * (nit_RCI2.GetPixel(k) - 127);

                    } // end loop over neighborhood

                    /** Set the result in the filtered part. */
                    itPSpf2[i].Set(tmp[i]);

                } // end loop over dimension i

                /** Increase all iterators. */
                ++nit_RCI;
                ++nit_RCI2;
                ++nit_RCI1;
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itPSpf2[i];
                    for (unsigned int j = 0; j < ImageDimension; j++) {
                        ++nitPSp2[i][j];
                    }
                }
            } // end while
        }   // end if do PrincipleStrainPart2

        nit_RCI.GoToBegin();
        nit_RCI1.GoToBegin();
        nit_RCI2.GoToBegin();
        if (this->m_CalculatePrincipleStrainPart3) {
            while (!itPSpf3[0].IsAtEnd()) {
                /** Create and reset tmp with zeros. */
                std::vector<double> tmp(ImageDimension, 0.0);

                /** Loop over all dimensions. */
                for (unsigned int i = 0; i < ImageDimension; i++) {
                    /** Loop over the neighborhood. */
                    for (unsigned int k = 0; k < neighborhoodSize; ++k) {
                        /** Calculation of the inner product. */
                        tmp[i] += Operator_A.GetElement(k)      // FA *
                                  * nitPSp3[i][0].GetPixel(k)          // subpart[ i ][ 0 ]
                                  * nit_RCI.GetPixel(k)
                                  * (nit_RCI1.GetPixel(k) - 127)
                                  * (nit_RCI2.GetPixel(k) - 127)
                                  * (nit_RCI1.GetPixel(k) - 127)
                                  * (nit_RCI2.GetPixel(k) - 127);

                        tmp[i] += Operator_B.GetElement(k)      // FB *
                                  * nitPSp3[i][1].GetPixel(k)          // subpart[ i ][ 1 ]
                                  * nit_RCI.GetPixel(k)
                                  * (nit_RCI1.GetPixel(k) - 127)
                                  * (nit_RCI2.GetPixel(k) - 127)
                                  * (nit_RCI1.GetPixel(k) - 127)
                                  * (nit_RCI2.GetPixel(k) - 127);

                    } // end loop over neighborhood

                    /** Set the result in the filtered part. */
                    itPSpf2[i].Set(tmp[i]);

                } // end loop over dimension i

                /** Increase all iterators. */
                ++nit_RCI;
                ++nit_RCI1;
                ++nit_RCI2;


                for (unsigned int i = 0; i < ImageDimension; i++) {
                    ++itPSpf3[i];
                    for (unsigned int j = 0; j < ImageDimension; j++) {
                        ++nitPSp3[i][j];
                    }
                }
            } // end while
        }   // end if do PrincipleStrainPart3



        /** TASK 8:
   * Add it all to create the final derivative images.
   ************************************************************************* */

        /** Create derivative images, each holding a component of the vector field. */
        std::vector <CoefficientImagePointer> derivativeImages(ImageDimension);
        for (unsigned int i = 0; i < ImageDimension; i++) {
            derivativeImages[i] = CoefficientImageType::New();
            derivativeImages[i]->SetRegions(inputImages[i]->GetLargestPossibleRegion());
            derivativeImages[i]->Allocate();
        }

        /** Create iterators over the derivative images. */
        std::vector <CoefficientImageIteratorType> itDIs(ImageDimension);
        for (unsigned int i = 0; i < ImageDimension; i++) {
            itDIs[i] = CoefficientImageIteratorType(derivativeImages[i],
                                                    derivativeImages[i]->GetLargestPossibleRegion());
            itDIs[i].GoToBegin();

            itGSpf1[i].GoToBegin();
            itGSpf2[i].GoToBegin();
            itPSpf1[i].GoToBegin();
            itPSpf2[i].GoToBegin();
            itPSpf3[i].GoToBegin();
        }

        /** Do the addition. */
        // NOTE: unlike the values, for the derivatives weight * derivative is returned.

        MeasureType gradMagGS1 = NumericTraits<MeasureType>::Zero;
        MeasureType gradMagGS2 = NumericTraits<MeasureType>::Zero;

        MeasureType gradMagPS1 = NumericTraits<MeasureType>::Zero;
        MeasureType gradMagPS2 = NumericTraits<MeasureType>::Zero;
        MeasureType gradMagPS3 = NumericTraits<MeasureType>::Zero;

        double rigidityCoefficientSumSqr = rigidityCoefficientSum * rigidityCoefficientSum;
        while (!itDIs[0].IsAtEnd()) {
            for (unsigned int i = 0; i < ImageDimension; i++) {
                ScalarType tmpDIs = NumericTraits<ScalarType>::Zero;

                /** Compute gradient magnitude of OC. */
                ScalarType tmpGS1 = this->m_GeneralStrainPartWeight * m_horizontalDirectionWeight * itGSpf1[i].Get();
                gradMagGS1 += tmpGS1 * tmpGS1 / rigidityCoefficientSumSqr;

                ScalarType tmpGS2 = this->m_GeneralStrainPartWeight * m_verticalDirectionWeight * itGSpf2[i].Get();
                gradMagGS2 += tmpGS2 * tmpGS2 / rigidityCoefficientSumSqr;


                /** Compute gradient magnitude of PC. */
                ScalarType tmpPS1 =
                        this->m_PrincipleStrainPartWeight * m_horizontalDirectionWeight * m_horizontalDirectionWeight *
                        itPSpf1[i].Get();
                gradMagPS1 += tmpPS1 * tmpPS1 / rigidityCoefficientSumSqr;

                ScalarType tmpPS2 =
                        this->m_PrincipleStrainPartWeight * m_verticalDirectionWeight * m_verticalDirectionWeight *
                        itPSpf2[i].Get();
                gradMagPS2 += tmpPS2 * tmpPS2 / rigidityCoefficientSumSqr;

                ScalarType tmpPS3 =
                        this->m_PrincipleStrainPartWeight * m_horizontalDirectionWeight * m_verticalDirectionWeight *
                        itPSpf3[i].Get();
                gradMagPS3 += tmpPS3 * tmpPS3 / rigidityCoefficientSumSqr;

                /** Compute derivative contribution. */

                if (this->m_UseGeneralStrainPart1) {
                    tmpDIs += tmpGS1;
                }
                if (this->m_UseGeneralStrainPart2) {
                    tmpDIs += tmpGS2;
                }
                if (this->m_UsePrincipleStrainPart1) {
                    tmpDIs += tmpPS1;
                }
                if (this->m_UsePrincipleStrainPart2) {
                    tmpDIs += tmpPS2;
                }
                if (this->m_UsePrincipleStrainPart3) {
                    tmpDIs += tmpPS3;
                }
                itDIs[i].Set(tmpDIs);

                /** Update iterators. */
                ++itDIs[i];
                ++itGSpf1[i];
                ++itGSpf2[i];
                ++itPSpf1[i];
                ++itPSpf2[i];
                ++itPSpf3[i];
            }
        } // end while

        /** Set the gradient magnitudes of the several terms. */

        this->m_GeneralStrainPartGradientMagnitude1 = vcl_sqrt(gradMagGS1);
        this->m_GeneralStrainPartGradientMagnitude2 = vcl_sqrt(gradMagGS2);

        this->m_PrincipleStrainPartGradientMagnitude1 = vcl_sqrt(gradMagPS1);
        this->m_PrincipleStrainPartGradientMagnitude2 = vcl_sqrt(gradMagPS2);
        this->m_PrincipleStrainPartGradientMagnitude3 = vcl_sqrt(gradMagPS3);

        /** Rearrange to create a derivative. */
        unsigned int j = 0;
        for (unsigned int i = 0; i < ImageDimension; i++) {
            itDIs[i].GoToBegin();
            while (!itDIs[i].IsAtEnd()) {
                derivative[j] = itDIs[i].Get() / rigidityCoefficientSum;
                ++itDIs[i];
                j++;
            } // end while
        }   // end for

        //std::cout << value << std::endl;


    } // end GetValueAndDerivative()


/**
 * ********************* PrintSelf ******************************
 */

    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::PrintSelf(std::ostream &os, Indent indent) const {
        /** Call the superclass' PrintSelf. */
        Superclass::PrintSelf(os, indent);

        /** Add debugging information. */


        os << indent << "GeneralStrainPartWeight: "
           << this->m_GeneralStrainPartWeight << std::endl;

        os << indent << "PrincipleStrainPartWeight: "
           << this->m_PrincipleStrainPartWeight << std::endl;

        os << indent << "horizontalDirectionWeight: "
           << this->m_horizontalDirectionWeight << std::endl;

        os << indent << "verticalDirectionWeight: "
           << this->m_verticalDirectionWeight << std::endl;

        os << indent << "RigidityCoefficientImage: "
           << this->m_RigidityCoefficientImage << std::endl;
        os << indent << "BSplineTransform: "
           << this->m_BSplineTransform << std::endl;
        os << indent << "m_StrainEnergyPenaltyTermValue : "
           << this->m_StrainEnergyPenaltyTermValue << std::endl;

        os << indent << "GeneralStrainPartValue1: "
           << this->m_GeneralStrainPartValue1 << std::endl;
        os << indent << "GeneralStrainPartValue2: "
           << this->m_GeneralStrainPartValue2 << std::endl;

        os << indent << "PrincipleStrainPartValue1: "
           << this->m_PrincipleStrainPartValue1 << std::endl;
        os << indent << "PrincipleStrainPartValue2: "
           << this->m_PrincipleStrainPartValue2 << std::endl;
        os << indent << "PrincipleStrainPartValue3: "
           << this->m_PrincipleStrainPartValue3 << std::endl;

        os << indent << "GeneralStrainPartGradientMagnitude1: "
           << this->m_GeneralStrainPartGradientMagnitude1 << std::endl;
        os << indent << "GeneralStrainPartGradientMagnitude2: "
           << this->m_GeneralStrainPartGradientMagnitude2 << std::endl;
        os << indent << "PrincipleStrainPartGradientMagnitude1: "
           << this->m_PrincipleStrainPartGradientMagnitude1 << std::endl;
        os << indent << "PrincipleStrainPartGradientMagnitude2: "
           << this->m_PrincipleStrainPartGradientMagnitude2 << std::endl;
        os << indent << "PrincipleStrainPartGradientMagnitude3: "
           << this->m_PrincipleStrainPartGradientMagnitude3 << std::endl;

        os << indent << "UseGeneralStrainPart1: "
           << this->m_UseGeneralStrainPart1 << std::endl;
        os << indent << "UseGeneralStrainPart2: "
           << this->m_UseGeneralStrainPart2 << std::endl;

        os << indent << "UsePrincipleStrainPart1: "
           << this->m_UsePrincipleStrainPart1 << std::endl;
        os << indent << "UsePrincipleStrainPart2: "
           << this->m_UsePrincipleStrainPart2 << std::endl;
        os << indent << "UsePrincipleStrainPart3: "
           << this->m_UsePrincipleStrainPart3 << std::endl;

        os << indent << "CalculateGeneralStrainPart1: "
           << this->m_CalculateGeneralStrainPart1 << std::endl;
        os << indent << "CalculateGeneralStrainPart2: "
           << this->m_CalculateGeneralStrainPart2 << std::endl;

        os << indent << "CalculatePrincipleStrainPart1: "
           << this->m_CalculatePrincipleStrainPart1 << std::endl;
        os << indent << "CalculatePrincipleStrainPart2: "
           << this->m_CalculatePrincipleStrainPart2 << std::endl;
        os << indent << "CalculatePrincipleStrainPart3: "
           << this->m_CalculatePrincipleStrainPart3 << std::endl;


    } // end PrintSelf()


/**
 * ************************ Create1DOperator *********************
 */
    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::Create1DOperator(
            NeighborhoodType &F,
            const std::string &WhichF,
            const unsigned int WhichDimension,
            const CoefficientImageSpacingType &spacing) const {
        /** Create an operator size and set it in the operator. */
        NeighborhoodSizeType r;
        r.Fill(NumericTraits<unsigned int>::ZeroValue());
        r[WhichDimension - 1] = 1;
        F.SetRadius(r);

        /** Get the image spacing factors that we are going to use. */
        std::vector<double> s(ImageDimension);
        for (unsigned int i = 0; i < ImageDimension; i++) {
            s[i] = spacing[i];
        }

        /** Create the required operator (neighborhood), depending on
   * WhichF. The operator is either 3x1 or 1x3 in 2D and
   * either 3x1x1 or 1x3x1 or 1x1x3 in 3D.
   */
        if (WhichF == "FA_xi" && WhichDimension == 1) {

            F[0] = -0.5 / s[0];
            F[1] = 0.0;
            F[2] = 0.5 / s[0];
        } else if (WhichF == "FA_xi" && WhichDimension == 2) {

            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FA_xi" && WhichDimension == 3) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FB_xi" && WhichDimension == 1) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FB_xi" && WhichDimension == 2) {
            F[0] = -0.5 / s[1];
            F[1] = 0.0;
            F[2] = 0.5 / s[1];
        } else if (WhichF == "FB_xi" && WhichDimension == 3) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FC_xi" && WhichDimension == 1) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FC_xi" && WhichDimension == 2) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FC_xi" && WhichDimension == 3) {
            F[0] = -0.5 / s[2];
            F[1] = 0.0;
            F[2] = 0.5 / s[2];
        } else if (WhichF == "FD_xi" && WhichDimension == 1) {

            F[0] = 0.5 / (s[0] * s[0]);
            F[1] = -1.0 / (s[0] * s[0]);
            F[2] = 0.5 / (s[0] * s[0]);
        } else if (WhichF == "FD_xi" && WhichDimension == 2) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FD_xi" && WhichDimension == 3) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FE_xi" && WhichDimension == 1) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FE_xi" && WhichDimension == 2) {
            F[0] = 0.5 / (s[1] * s[1]);
            F[1] = -1.0 / (s[1] * s[1]);
            F[2] = 0.5 / (s[1] * s[1]);
        } else if (WhichF == "FE_xi" && WhichDimension == 3) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FF_xi" && WhichDimension == 1) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FF_xi" && WhichDimension == 2) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FF_xi" && WhichDimension == 3) {
            F[0] = 0.5 / (s[2] * s[2]);
            F[1] = -1.0 / (s[2] * s[2]);
            F[2] = 0.5 / (s[2] * s[2]);
        } else if (WhichF == "FG_xi" && WhichDimension == 1) {
            F[0] = -0.5 / (s[0] * s[1]);
            F[1] = 0.0;
            F[2] = 0.5 / (s[0] * s[1]);
        } else if (WhichF == "FG_xi" && WhichDimension == 2) {
            F[0] = -0.5 / (s[0] * s[1]);
            F[1] = 0.0;
            F[2] = 0.5 / (s[0] * s[1]);
        } else if (WhichF == "FG_xi" && WhichDimension == 3) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FH_xi" && WhichDimension == 1) {
            F[0] = -0.5 / (s[0] * s[2]);
            F[1] = 0.0;
            F[2] = 0.5 / (s[0] * s[2]);
        } else if (WhichF == "FH_xi" && WhichDimension == 2) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FH_xi" && WhichDimension == 3) {
            F[0] = -0.5 / (s[0] * s[2]);
            F[1] = 0.0;
            F[2] = 0.5 / (s[0] * s[2]);
        } else if (WhichF == "FI_xi" && WhichDimension == 1) {
            F[0] = 1.0 / 6.0;
            F[1] = 4.0 / 6.0;
            F[2] = 1.0 / 6.0;
        } else if (WhichF == "FI_xi" && WhichDimension == 2) {
            F[0] = -0.5 / (s[1] * s[2]);
            F[1] = 0.0;
            F[2] = 0.5 / (s[1] * s[2]);
        } else if (WhichF == "FI_xi" && WhichDimension == 3) {
            F[0] = -0.5 / (s[1] * s[2]);
            F[1] = 0.0;
            F[2] = 0.5 / (s[1] * s[2]);
        } else {
            /** Throw an exception. */
            itkExceptionMacro( << "Can not create this type of operator." );
        }

    } // end Create1DOperator()

/////////////////////////////////////////////////////////
/**
 * ************************** FilterSeparable ********************
 */

    template<class TFixedImage, class TScalarType>
    typename TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>::CoefficientImagePointer
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::FilterSeparable(
            const CoefficientImageType *image,
            const std::vector <NeighborhoodType> &Operators) const {
        /** Create filters, supply them with boundary conditions and operators. */
        std::vector<typename NOIFType::Pointer> filters(ImageDimension);
        for (unsigned int i = 0; i < ImageDimension; i++) {
            filters[i] = NOIFType::New();
            filters[i]->SetOperator(Operators[i]);
        }

        /** Set up the mini-pipline. */
        filters[0]->SetInput(image);
        for (unsigned int i = 1; i < ImageDimension; i++) {
            filters[i]->SetInput(filters[i - 1]->GetOutput());
        }

        /** Execute the mini-pipeline. */
        filters[ImageDimension - 1]->Update();

        /** Return the filtered image. */
        return filters[ImageDimension - 1]->GetOutput();

    } // end FilterSeparable()


/**
 * ************************ CreateNDOperator *********************
 */

    template<class TFixedImage, class TScalarType>
    void
    TransformDirectionalStrainEnergyPenaltyTerm<TFixedImage, TScalarType>
    ::CreateNDOperator(
            NeighborhoodType &F,
            const std::string &WhichF,
            const CoefficientImageSpacingType &spacing) const {
        /** Create an operator size and set it in the operator. */
        NeighborhoodSizeType r;
        r.Fill(1);
        F.SetRadius(r);

        /** Get the image spacing factors that we are going to use. */
        std::vector<double> s(ImageDimension);
        for (unsigned int i = 0; i < ImageDimension; i++) {
            s[i] = spacing[i];
        }

        /** Create the required operator (neighborhood), depending on
   * WhichF. The operator is either 3x3 in 2D or 3x3x3 in 3D.
   */
        if (WhichF == "FA") {
            if (ImageDimension == 2) {
                F[0] = 1.0 / 12.0 / s[0];
                F[1] = 0.0;
                F[2] = -1.0 / 12.0 / s[0];
                F[3] = 1.0 / 3.0 / s[0];
                F[4] = 0.0;
                F[5] = -1.0 / 3.0 / s[0];
                F[6] = 1.0 / 12.0 / s[0];
                F[7] = 0.0;
                F[8] = -1.0 / 12.0 / s[0];
            }

        } else if (WhichF == "FB") {
            if (ImageDimension == 2) {
                F[0] = 1.0 / 12.0 / s[1];
                F[1] = 1.0 / 3.0 / s[1];
                F[2] = 1.0 / 12.0 / s[1];
                F[3] = 0.0;
                F[4] = 0.0;
                F[5] = 0.0;
                F[6] = -1.0 / 12.0 / s[1];
                F[7] = -1.0 / 3.0 / s[1];
                F[8] = -1.0 / 12.0 / s[1];
            }

        } else if (WhichF == "FD") {
            if (ImageDimension == 2) {
                double sp = s[0] * s[0];
                F[0] = 1.0 / 12.0 / sp;
                F[1] = -1.0 / 6.0 / sp;
                F[2] = 1.0 / 12.0 / sp;
                F[3] = 1.0 / 3.0 / sp;
                F[4] = -2.0 / 3.0 / sp;
                F[5] = 1.0 / 3.0 / sp;
                F[6] = 1.0 / 12.0 / sp;
                F[7] = -1.0 / 6.0 / sp;
                F[8] = 1.0 / 12.0 / sp;
            }

        } else if (WhichF == "FE") {
            if (ImageDimension == 2) {
                double sp = s[1] * s[1];
                F[0] = 1.0 / 12.0 / sp;
                F[1] = 1.0 / 3.0 / sp;
                F[2] = 1.0 / 12.0 / sp;
                F[3] = -1.0 / 6.0 / sp;
                F[4] = -2.0 / 3.0 / sp;
                F[5] = -1.0 / 6.0 / sp;
                F[6] = 1.0 / 12.0 / sp;
                F[7] = 1.0 / 3.0 / sp;
                F[8] = 1.0 / 12.0 / sp;
            }

        } else if (WhichF == "FG") {
            if (ImageDimension == 2) {
                double sp = s[0] * s[1];
                F[0] = 1.0 / 4.0 / sp;
                F[1] = 0.0;
                F[2] = -1.0 / 4.0 / sp;
                F[3] = 0.0;
                F[4] = 0.0;
                F[5] = 0.0;
                F[6] = -1.0 / 4.0 / sp;
                F[7] = 0.0;
                F[8] = 1.0 / 4.0 / sp;
            }

        } else {
            /** Throw an exception. */
            itkExceptionMacro( << "Can not create this type of operator." );
        }

    } // end CreateNDOperator()


} // end namespace itk

#endif // #ifndef __itkTransformDirectionalStrainEnergyPenaltyTerm_hxx



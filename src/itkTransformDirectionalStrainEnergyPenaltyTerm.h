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

#ifndef __itkTransformDirectionalStrainEnergyPenaltyTerm_h
#define __itkTransformDirectionalStrainEnergyPenaltyTerm_h

#include "itkTransformPenaltyTerm.h"
#include "itkImageGridSampler.h"

/** Needed for the check of a B-spline transform. */
#include "itkAdvancedBSplineDeformableTransform.h"
#include "itkAdvancedCombinationTransform.h"

/** Needed for the filtering of the B-spline coefficients. */
#include "itkNeighborhood.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodOperatorImageFilter.h"
#include "itkNeighborhoodIterator.h"

/** Include stuff needed for the construction of the coefficient images. */
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkImageRegionIterator.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

namespace itk {

    template<class TFixedImage, class TScalarType>
    class TransformDirectionalStrainEnergyPenaltyTerm :
            public TransformPenaltyTerm<TFixedImage, TScalarType> {
    public:

        /** Standard ITK stuff. */
        typedef TransformDirectionalStrainEnergyPenaltyTerm Self;
        typedef TransformPenaltyTerm <
        TFixedImage, TScalarType> Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro( Self );

        /** Run-time type information (and related methods). */
        itkTypeMacro( TransformDirectionalStrainEnergyPenaltyTerm, TransformPenaltyTerm
        );

        /** Typedefs inherited from the superclass. */
        typedef typename Superclass::CoordinateRepresentationType CoordinateRepresentationType;
        typedef typename Superclass::MovingImageType MovingImageType;
        typedef typename Superclass::MovingImagePixelType MovingImagePixelType;
        typedef typename Superclass::MovingImagePointer MovingImagePointer;
        typedef typename Superclass::MovingImageConstPointer MovingImageConstPointer;
        typedef typename Superclass::FixedImageType FixedImageType;
        typedef typename Superclass::FixedImagePointer FixedImagePointer;
        typedef typename Superclass::FixedImageConstPointer FixedImageConstPointer;
        typedef typename Superclass::FixedImageRegionType FixedImageRegionType;
        typedef typename Superclass::TransformType TransformType;
        typedef typename Superclass::TransformPointer TransformPointer;
        typedef typename Superclass::InputPointType InputPointType;
        typedef typename Superclass::OutputPointType OutputPointType;
        typedef typename Superclass::TransformParametersType TransformParametersType;
        typedef typename Superclass::TransformJacobianType TransformJacobianType;
        typedef typename Superclass::NumberOfParametersType NumberOfParametersType;
        typedef typename Superclass::InterpolatorType InterpolatorType;
        typedef typename Superclass::InterpolatorPointer InterpolatorPointer;
        typedef typename Superclass::RealType RealType;
        typedef typename Superclass::GradientPixelType GradientPixelType;
        typedef typename Superclass::GradientImageType GradientImageType;
        typedef typename Superclass::GradientImagePointer GradientImagePointer;
        typedef typename Superclass::GradientImageFilterType GradientImageFilterType;
        typedef typename Superclass::GradientImageFilterPointer GradientImageFilterPointer;
        typedef typename Superclass::FixedImageMaskType FixedImageMaskType;
        typedef typename Superclass::FixedImageMaskPointer FixedImageMaskPointer;
        typedef typename Superclass::MovingImageMaskType MovingImageMaskType;
        typedef typename Superclass::MovingImageMaskPointer MovingImageMaskPointer;
        typedef typename Superclass::MeasureType MeasureType;
        typedef typename Superclass::DerivativeType DerivativeType;
        typedef typename Superclass::DerivativeValueType DerivativeValueType;
        typedef typename Superclass::ParametersType ParametersType;
        typedef typename Superclass::FixedImagePixelType FixedImagePixelType;
        typedef typename Superclass::ImageSampleContainerType ImageSampleContainerType;
        typedef typename Superclass::ImageSampleContainerPointer ImageSampleContainerPointer;
        typedef typename Superclass::ScalarType ScalarType;
        typedef typename Superclass::ThreaderType ThreaderType;
        typedef typename Superclass::ThreadInfoType ThreadInfoType;

        /** Typedefs from the AdvancedTransform. */
        typedef typename Superclass::SpatialJacobianType SpatialJacobianType;
        typedef typename Superclass
        ::JacobianOfSpatialJacobianType JacobianOfSpatialJacobianType;
        typedef typename Superclass::SpatialHessianType SpatialHessianType;
        typedef typename Superclass
        ::JacobianOfSpatialHessianType JacobianOfSpatialHessianType;
        typedef typename Superclass::InternalMatrixType InternalMatrixType;
        typedef typename Superclass::HessianValueType HessianValueType;
        typedef typename Superclass::HessianType HessianType;

        /** Define the dimension. */
        itkStaticConstMacro( FixedImageDimension,
        unsigned int, FixedImageType::ImageDimension );
        itkStaticConstMacro( MovingImageDimension,
        unsigned int, FixedImageType::ImageDimension );
        itkStaticConstMacro( ImageDimension,
        unsigned int, FixedImageType::ImageDimension );

        /** Initialize the penalty term. */
        virtual void Initialize(void) throw(ExceptionObject);

        /** Typedef's for B-spline transform. */
        typedef AdvancedBSplineDeformableTransform<ScalarType,
                FixedImageDimension, 3> BSplineTransformType;
        typedef typename BSplineTransformType::Pointer BSplineTransformPointer;
        typedef typename BSplineTransformType::SpacingType GridSpacingType;
        typedef typename BSplineTransformType::ImageType CoefficientImageType;
        typedef typename CoefficientImageType::Pointer CoefficientImagePointer;
        typedef typename CoefficientImageType::SpacingType CoefficientImageSpacingType;
        typedef typename BSplineTransformType::ImageType HorizontalVectorImageType;
        typedef typename HorizontalVectorImageType::Pointer HorizontalVectorImagePointer;
        typedef typename HorizontalVectorImageType::SpacingType HorizontalVectorImageSpacingType;
        typedef typename BSplineTransformType::ImageType VerticalVectorImageType;
        typedef typename VerticalVectorImageType::Pointer VerticalVectorImagePointer;
        typedef typename VerticalVectorImageType::SpacingType VerticalVectorImageSpacingType;
        typedef AdvancedCombinationTransform <ScalarType,
        FixedImageDimension> CombinationTransformType;

        /** Typedef support for neighborhoods, filters, etc. */
        typedef Neighborhood<ScalarType,
                itkGetStaticConstMacro(FixedImageDimension)> NeighborhoodType;

        typedef typename NeighborhoodType::SizeType NeighborhoodSizeType;

        typedef ImageRegionIterator <CoefficientImageType> CoefficientImageIteratorType;
        typedef NeighborhoodOperatorImageFilter <
        CoefficientImageType, CoefficientImageType> NOIFType;
        typedef NeighborhoodIterator <CoefficientImageType> NeighborhoodIteratorType;
        typedef typename NeighborhoodIteratorType::RadiusType RadiusType;


        /***~^ ^~****/
        typedef ImageRegionIterator <HorizontalVectorImageType> HorizontalVectorImageIteratorType;
        typedef NeighborhoodOperatorImageFilter <
        HorizontalVectorImageType, HorizontalVectorImageType> NOIFType1;
        typedef NeighborhoodIterator <HorizontalVectorImageType> NeighborhoodIteratorType1;
        typedef typename NeighborhoodIteratorType::RadiusType RadiusType1;
        /******/
        typedef ImageRegionIterator <VerticalVectorImageType> VerticalVectorImageIteratorType;
        typedef NeighborhoodOperatorImageFilter <
        VerticalVectorImageType, VerticalVectorImageType> NOIFType2;
        typedef NeighborhoodIterator <VerticalVectorImageType> NeighborhoodIteratorType2;
        typedef typename NeighborhoodIteratorType::RadiusType RadiusType2;
        /***~^ ^~****/


        /** Typedef's for the construction of the coefficient images. */
        typedef CoefficientImageType RigidityImageType;

        typedef typename RigidityImageType::Pointer RigidityImagePointer;
        typedef typename RigidityImageType::PixelType RigidityPixelType;
        typedef typename RigidityImageType::RegionType RigidityImageRegionType;
        typedef typename RigidityImageType::IndexType RigidityImageIndexType;
        typedef typename RigidityImageType::PointType RigidityImagePointType;

        typedef ImageRegionIterator <RigidityImageType> RigidityImageIteratorType;
        typedef BinaryBallStructuringElement<
                RigidityPixelType,
                itkGetStaticConstMacro(FixedImageDimension)> StructuringElementType;
        typedef typename StructuringElementType::RadiusType SERadiusType;
        typedef GrayscaleDilateImageFilter <
        RigidityImageType, RigidityImageType,
        StructuringElementType> DilateFilterType;
        typedef typename DilateFilterType::Pointer DilateFilterPointer;
        typedef HorizontalVectorImageType VectorImageType1;
        typedef VerticalVectorImageType VectorImageType2;

        typedef typename VectorImageType1::Pointer VectorImagePointer1;
        typedef typename VectorImageType1::PixelType VectorPixelType1;
        typedef typename VectorImageType1::RegionType VectorImageRegionType1;
        typedef typename VectorImageType1::IndexType VectorImageIndexType1;
        typedef typename VectorImageType1::PointType VectorImagePointType1;
        typedef ImageRegionIterator <VectorImageType1> VectorImageIteratorType1;
        typedef BinaryBallStructuringElement<
                VectorPixelType1,
                itkGetStaticConstMacro(FixedImageDimension)> StructuringElementType1;
        typedef GrayscaleDilateImageFilter <
        VectorImageType1, VectorImageType1,
        StructuringElementType1> DilateFilterType1;


        typedef NearestNeighborInterpolateImageFunction<
                HorizontalVectorImageType, double> NearestNeighborInterpolatorType1;

        typedef ResampleImageFilter <HorizontalVectorImageType, HorizontalVectorImageType>
                ResampleFilterType;

        typedef typename VectorImageType2::Pointer VectorImagePointer2;
        typedef typename VectorImageType2::PixelType VectorPixelType2;
        typedef typename VectorImageType2::RegionType VectorImageRegionType2;
        typedef typename VectorImageType2::IndexType VectorImageIndexType2;
        typedef typename VectorImageType2::PointType VectorImagePointType2;
        typedef ImageRegionIterator <VectorImageType2> VectorImageIteratorType2;
        typedef BinaryBallStructuringElement<
                VectorPixelType2,
                itkGetStaticConstMacro(FixedImageDimension)> StructuringElementType2;
        typedef GrayscaleDilateImageFilter <
        VectorImageType2, VectorImageType2,
        StructuringElementType2> DilateFilterType2;

        typedef NearestNeighborInterpolateImageFunction<
                VectorImageType2, double> NearestNeighborInterpolatorType2;
        typedef typename NearestNeighborInterpolatorType2::Pointer nearestNeighborInterpolatorType2;

        /** Check stuff. */
        void CheckUseAndCalculationBooleans(void);

        /** The GetValue()-method returns the rigid penalty value. */
        virtual MeasureType GetValue(
                const ParametersType &parameters) const;

        /** The GetDerivative()-method returns the rigid penalty derivative. */
        virtual void GetDerivative(
                const ParametersType &parameters,
                DerivativeType &derivative) const;

        /** Contains calls from GetValueAndDerivative that are thread-unsafe. */
        virtual void BeforeThreadedGetValueAndDerivative(
                const TransformParametersType &parameters) const;

        /** The GetValueAndDerivative()-method returns the rigid penalty value and its derivative. */
        virtual void GetValueAndDerivative(
                const ParametersType &parameters,
                MeasureType &value,
                DerivativeType &derivative) const;

        /** Set the B-spline transform in this class.
         * This class expects a BSplineTransform! It is not suited for others.
         */
        itkSetObjectMacro( BSplineTransform, BSplineTransformType
        );

        itkSetClampMacro( PreSetParameter, ScalarType,
        0.0,

        NumericTraits<ScalarType>::max()

        );
        itkGetMacro( PreSetParameter, ScalarType
        );

        itkSetClampMacro( GeneralStrainPartWeight, ScalarType,
        0.0,

        NumericTraits<ScalarType>::max()

        );
        itkGetMacro( GeneralStrainPartWeight, ScalarType
        );

        /** Set/Get the weight of the PrincipleStrain condition part. */
        itkSetClampMacro( PrincipleStrainPartWeight, ScalarType,
        0.0,

        NumericTraits<ScalarType>::max()

        );
        itkGetMacro( PrincipleStrainPartWeight, ScalarType
        );


        itkSetClampMacro( horizontalDirectionWeight, ScalarType,
        0.0,

        NumericTraits<ScalarType>::max()

        );
        itkGetMacro( horizontalDirectionWeight, ScalarType
        );

        itkSetClampMacro( verticalDirectionWeight, ScalarType,
        0.0,

        NumericTraits<ScalarType>::max()

        );
        itkGetMacro( verticalDirectionWeight, ScalarType
        );

        /** Set the usage of the GeneralStrainPart. */
        itkSetMacro( UseGeneralStrainPart1,
        bool );
        itkSetMacro( UseGeneralStrainPart2,
        bool );

        /** Set the usage of the PrincipleStrain condition part. */
        itkSetMacro( UsePrincipleStrainPart1,
        bool );
        itkSetMacro( UsePrincipleStrainPart2,
        bool );
        itkSetMacro( UsePrincipleStrainPart3,
        bool );

        /** Set the calculation of the GeneralStrainPart,
         * even if we don't use it.
         */
        itkSetMacro( CalculateGeneralStrainPart1,
        bool );
        itkSetMacro( CalculateGeneralStrainPart2,
        bool );
        /** Set the calculation of the properness condition part.,
         * even if we don't use it.
         */
        itkSetMacro( CalculatePrincipleStrainPart1,
        bool );
        itkSetMacro( CalculatePrincipleStrainPart2,
        bool );
        itkSetMacro( CalculatePrincipleStrainPart3,
        bool );

        /** Get the value of the GeneralStrainPart. */
        itkGetConstReferenceMacro( GeneralStrainPartValue1, MeasureType
        );
        itkGetConstReferenceMacro( GeneralStrainPartValue2, MeasureType
        );
        /** Get the value of the properness condition. */
        itkGetConstReferenceMacro( PrincipleStrainPartValue1, MeasureType
        );
        itkGetConstReferenceMacro( PrincipleStrainPartValue2, MeasureType
        );
        itkGetConstReferenceMacro( PrincipleStrainPartValue3, MeasureType
        );

        /** Get the gradient magnitude of the GeneralStrainPart. */
        itkGetConstReferenceMacro( GeneralStrainPartGradientMagnitude1, MeasureType
        );
        itkGetConstReferenceMacro( GeneralStrainPartGradientMagnitude2, MeasureType
        );
        /** Get the gradient magnitude of the PrincipleStrain condition. */
        itkGetConstReferenceMacro( PrincipleStrainPartGradientMagnitude1, MeasureType
        );
        itkGetConstReferenceMacro( PrincipleStrainPartGradientMagnitude2, MeasureType
        );
        itkGetConstReferenceMacro( PrincipleStrainPartGradientMagnitude3, MeasureType
        );
        /** Get the value of the total StrainEnergy penalty term. */
        //itkGetConstReferenceMacro( RigidityPenaltyTermValue, MeasureType );

        /** Set if the RigidityImage's are dilated. */
        itkSetMacro( DilateRigidityImages,
        bool );
        itkSetMacro( DilateVectorImages,
        bool );

        /** Set the DilationRadiusMultiplier. */
        itkSetClampMacro( DilationRadiusMultiplier, CoordinateRepresentationType,
        0.1,

        NumericTraits<CoordinateRepresentationType>::max()

        );

        /** Set the fixed coefficient image. */
        itkSetObjectMacro( FixedRigidityImage, RigidityImageType
        );
        itkSetObjectMacro( HorizontalDirectionVectorImage, VectorImageType1
        );

        /** Set the moving coefficient image. */
        itkSetObjectMacro( MovingRigidityImage, RigidityImageType
        );

        itkSetObjectMacro( VerticalDirectionVectorImage, VectorImageType2
        );

        /** Set to use the FixedRigidityImage or not. */
        itkSetMacro( UseFixedRigidityImage,
        bool );

        itkSetMacro( UseHorizontalDirectionVectorImage,
        bool );

        /** Set to use the MovingRigidityImage or not. */
        itkSetMacro( UseMovingRigidityImage,
        bool );

        itkSetMacro( UseVerticalDirectionVectorImage,
        bool );

        /** Function to fill the RigidityCoefficientImage every iteration. */
        void FillRigidityCoefficientImage(const ParametersType &parameters) const;

        void FillRigidityHorizontalVectorImage(const ParametersType &parameters) const;

        void FillRigidityVerticalVectorImage(const ParametersType &parameters) const;

    protected:

        /** The constructor. */
        TransformDirectionalStrainEnergyPenaltyTerm();

        /** The destructor. */
        virtual ~TransformDirectionalStrainEnergyPenaltyTerm() {}

        /** PrintSelf. */
        void PrintSelf(std::ostream &os, Indent indent) const;

    private:

        /** The private constructor. */
        TransformDirectionalStrainEnergyPenaltyTerm(const Self &); // purposely not implemented
        /** The private copy constructor. */
        void operator=(const Self &);            // purposely not implemented

        /** Internal function to dilate the rigidity images. */
        virtual void DilateRigidityImages(void);

        virtual void DilateVectorImages1(void);

        virtual void DilateVectorImages2(void);

        /** Private function used for the filtering. It creates 1D separable operators F. */
        void Create1DOperator(NeighborhoodType &F, const std::string &whichF,
                              const unsigned int WhichDimension, const CoefficientImageSpacingType &spacing) const;

        void Create1DOperator1(NeighborhoodType &F1, const std::string &whichF1,
                               const unsigned int WhichDimension1,
                               const HorizontalVectorImageSpacingType &spacing1) const;

        void Create1DOperator2(NeighborhoodType &F2, const std::string &whichF2,
                               const unsigned int WhichDimension2,
                               const VerticalVectorImageSpacingType &spacing2) const;

        /** Private function used for the filtering. It creates ND inseparable operators F. */
        void CreateNDOperator(NeighborhoodType &F, const std::string &whichF,
                              const CoefficientImageSpacingType &spacing) const;

        /** Private function used for the filtering. It performs 1D separable filtering. */
        CoefficientImagePointer FilterSeparable(const CoefficientImageType *,
                                                const std::vector <NeighborhoodType> &Operators) const;

        void CreateNDOperatorN1(NeighborhoodType &FN1, const std::string &whichFN1,
                                const HorizontalVectorImageSpacingType &spacingN1) const;

        HorizontalVectorImagePointer FilterSeparable1(const HorizontalVectorImageType *,
                                                      const std::vector <NeighborhoodType> &Operator1) const;

        void CreateNDOperatorN2(NeighborhoodType &FN2, const std::string &whichN2,
                                const VerticalVectorImageSpacingType &spacingN2) const;

        VerticalVectorImagePointer FilterSeparable2(const VerticalVectorImageType *,
                                                    const std::vector <NeighborhoodType> &Operators2) const;

        /** Member variables. */
        BSplineTransformPointer m_BSplineTransform;

        ScalarType m_PreSetParameter;
        ScalarType m_GeneralStrainPartWeight;
        ScalarType m_PrincipleStrainPartWeight;
        ScalarType m_horizontalDirectionWeight;
        ScalarType m_verticalDirectionWeight;

        mutable MeasureType m_StrainEnergyPenaltyTermValue;

        mutable MeasureType m_GeneralStrainPartValue1;
        mutable MeasureType m_GeneralStrainPartValue2;

        mutable MeasureType m_PrincipleStrainPartValue1;
        mutable MeasureType m_PrincipleStrainPartValue2;
        mutable MeasureType m_PrincipleStrainPartValue3;

        mutable MeasureType m_GeneralStrainPartGradientMagnitude1;
        mutable MeasureType m_GeneralStrainPartGradientMagnitude2;

        mutable MeasureType m_PrincipleStrainPartGradientMagnitude1;
        mutable MeasureType m_PrincipleStrainPartGradientMagnitude2;
        mutable MeasureType m_PrincipleStrainPartGradientMagnitude3;


        bool m_UseGeneralStrainPart1;
        bool m_UseGeneralStrainPart2;

        bool m_UsePrincipleStrainPart1;
        bool m_UsePrincipleStrainPart2;
        bool m_UsePrincipleStrainPart3;

        bool m_CalculateGeneralStrainPart1;
        bool m_CalculateGeneralStrainPart2;

        bool m_CalculatePrincipleStrainPart1;
        bool m_CalculatePrincipleStrainPart2;
        bool m_CalculatePrincipleStrainPart3;
        /** Rigidity image variables. */
        CoordinateRepresentationType m_DilationRadiusMultiplier;
        bool m_DilateRigidityImages;
        bool m_DilateVectorImages;

        mutable bool m_RigidityCoefficientImageIsFilled;
        mutable bool m_RigidityHorizontalVectorImageIsFilled;
        mutable bool m_RigidityVerticalVectorImageIsFilled;

        RigidityImagePointer m_FixedRigidityImage;
        RigidityImagePointer m_MovingRigidityImage;
        RigidityImagePointer m_RigidityCoefficientImage;

        VectorImagePointer1 m_HorizontalDirectionVectorImage;
        VectorImagePointer2 m_VerticalDirectionVectorImage;
        VectorImagePointer1 m_RigidityHorizontalVectorImage;
        VectorImagePointer2 m_RigidityVerticalVectorImage;

        std::vector <DilateFilterPointer> m_FixedRigidityImageDilation;
        std::vector <DilateFilterPointer> m_MovingRigidityImageDilation;
        RigidityImagePointer m_FixedRigidityImageDilated;
        RigidityImagePointer m_MovingRigidityImageDilated;
        bool m_UseFixedRigidityImage;
        bool m_UseMovingRigidityImage;

        std::vector <DilateFilterPointer> m_HorizontalDirectionVectorImageDilation;
        std::vector <DilateFilterPointer> m_VerticalDirectionVectorImageDilation;
        VectorImagePointer1 m_HorizontalDirectionVectorImageDilated;
        VectorImagePointer2 m_VerticalDirectionVectorImageDilated;
        bool m_UseHorizontalDirectionVectorImage;
        bool m_UseVerticalDirectionVectorImage;

    };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION

#include "itkTransformDirectionalStrainEnergyPenaltyTerm.hxx"

#endif

#endif // #ifndef __itkTransformDirectionalStrainEnergyPenaltyTerm_h

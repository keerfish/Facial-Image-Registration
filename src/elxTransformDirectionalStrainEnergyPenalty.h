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
#ifndef __elxTransformDirectionalStrainEnergyPenaltyTerm_H__
#define __elxTransformDirectionalStrainEnergyPenaltyTerm_H__

#include <DirectionalStrainEnergyPenalty/itkTransformDirectionalStrainEnergyPenaltyTerm.h>
#include "elxIncludes.h" // include first to avoid MSVS warning

namespace elastix {

#ifndef __elxTransformDirectionalStrainEnergyPenaltyTerm_H__
#define __elxTransformDirectionalStrainEnergyPenaltyTerm_H__

#include <DirectionalStrainEnergyPenalty/itkTransformDirectionalStrainEnergyPenaltyTerm.h>
#include "elxIncludes.h" // include first to avoid MSVS warning

    namespace elastix
    {

    template< class TElastix >
    class TransformDirectionalStrainEnergyPenalty :
      public
      itk::TransformDirectionalStrainEnergyPenaltyTerm<
      typename MetricBase< TElastix >::FixedImageType,
      double >,
      public MetricBase< TElastix >
    {
    public:

      /** Standard ITK-stuff. */
      typedef TransformDirectionalStrainEnergyPenalty Self;
      typedef itk::TransformDirectionalStrainEnergyPenaltyTerm<
      typename MetricBase< TElastix >::FixedImageType,
        double >                                            Superclass1;
      typedef MetricBase< TElastix >          Superclass2;
      typedef itk::SmartPointer< Self >       Pointer;
      typedef itk::SmartPointer< const Self > ConstPointer;

      /** Method for creation through the object factory. */
      itkNewMacro( Self );

      /** Run-time type information (and related methods). */
      itkTypeMacro( TransformDirectionalStrainEnergyPenalty, itk::TransformDirectionalStrainEnergyPenaltyTerm );

      /** Name of this class.
       * Use this name in the parameter file to select this specific metric. \n
       * example: <tt>(Metric "TransformDirectionalStrainEnergyPenalty")</tt>\n
       */
      elxClassNameMacro( "TransformDirectionalStrainEnergyPenalty" );

      /** Typedefs from the superclass. */
      typedef typename
        Superclass1::CoordinateRepresentationType CoordinateRepresentationType;
      typedef typename Superclass1::MovingImageType            MovingImageType;
      typedef typename Superclass1::MovingImagePixelType       MovingImagePixelType;
      typedef typename Superclass1::MovingImageConstPointer    MovingImageConstPointer;
      typedef typename Superclass1::FixedImageType             FixedImageType;
      typedef typename Superclass1::FixedImageConstPointer     FixedImageConstPointer;
      typedef typename Superclass1::FixedImageRegionType       FixedImageRegionType;
      typedef typename Superclass1::TransformType              TransformType;
      typedef typename Superclass1::TransformPointer           TransformPointer;
      typedef typename Superclass1::InputPointType             InputPointType;
      typedef typename Superclass1::OutputPointType            OutputPointType;
      typedef typename Superclass1::TransformParametersType    TransformParametersType;
      typedef typename Superclass1::TransformJacobianType      TransformJacobianType;
      typedef typename Superclass1::InterpolatorType           InterpolatorType;
      typedef typename Superclass1::InterpolatorPointer        InterpolatorPointer;
      typedef typename Superclass1::RealType                   RealType;
      typedef typename Superclass1::GradientPixelType          GradientPixelType;
      typedef typename Superclass1::GradientImageType          GradientImageType;
      typedef typename Superclass1::GradientImagePointer       GradientImagePointer;
      typedef typename Superclass1::GradientImageFilterType    GradientImageFilterType;
      typedef typename Superclass1::GradientImageFilterPointer GradientImageFilterPointer;
      typedef typename Superclass1::FixedImageMaskType         FixedImageMaskType;
      typedef typename Superclass1::FixedImageMaskPointer      FixedImageMaskPointer;
      typedef typename Superclass1::MovingImageMaskType        MovingImageMaskType;
      typedef typename Superclass1::MovingImageMaskPointer     MovingImageMaskPointer;
      typedef typename Superclass1::MeasureType                MeasureType;
      typedef typename Superclass1::DerivativeType             DerivativeType;
      typedef typename Superclass1::ParametersType             ParametersType;
      typedef typename Superclass1::FixedImagePixelType        FixedImagePixelType;
      typedef typename Superclass1::MovingImageRegionType      MovingImageRegionType;
      typedef typename Superclass1::ImageSamplerType           ImageSamplerType;
      typedef typename Superclass1::ImageSamplerPointer        ImageSamplerPointer;
      typedef typename Superclass1::ImageSampleContainerType   ImageSampleContainerType;
      typedef typename
        Superclass1::ImageSampleContainerPointer ImageSampleContainerPointer;
      typedef typename Superclass1::FixedImageLimiterType  FixedImageLimiterType;
      typedef typename Superclass1::MovingImageLimiterType MovingImageLimiterType;
      typedef typename
        Superclass1::FixedImageLimiterOutputType FixedImageLimiterOutputType;
      typedef typename
        Superclass1::MovingImageLimiterOutputType MovingImageLimiterOutputType;

      /** The fixed image dimension. */
      itkStaticConstMacro( FixedImageDimension, unsigned int,
        FixedImageType::ImageDimension );

      /** The moving image dimension. */
      itkStaticConstMacro( MovingImageDimension, unsigned int,
        MovingImageType::ImageDimension );

      /** Typedef's inherited from elastix. */
      typedef typename Superclass2::ElastixType          ElastixType;
      typedef typename Superclass2::ElastixPointer       ElastixPointer;
      typedef typename Superclass2::ConfigurationType    ConfigurationType;
      typedef typename Superclass2::ConfigurationPointer ConfigurationPointer;
      typedef typename Superclass2::RegistrationType     RegistrationType;
      typedef typename Superclass2::RegistrationPointer  RegistrationPointer;
      typedef typename Superclass2::ITKBaseType          ITKBaseType;

      /** Sets up a timer to measure the initialization time and
       * calls the Superclass' implementation.
       */
      virtual void Initialize( void ) throw ( itk::ExceptionObject );

      /**
       * Do some things before each resolution:
       * \li Set options for SelfHessian
       */
      virtual void BeforeEachResolution( void );

      virtual void BeforeRegistration( void );

        /**
         * Do some things after each iteration:
         * \li Print the OC, PC, LC parts of the rigidity term.
         */
        virtual void AfterEachIteration( void );

        /** This metric is advanced (so it has a sampling possibility), but it
         * purposely does not use samplers. The MetricBase class, however, issues
         * a warning if this is the case, so we overwrite that function.
         */
        virtual void SelectNewSamples( void ){}

    protected:

      /** The constructor. */
      TransformDirectionalStrainEnergyPenalty(){}

      /** The destructor. */
      virtual ~TransformDirectionalStrainEnergyPenalty() {}

    private:

      /** The private constructor. */
      TransformDirectionalStrainEnergyPenalty( const Self & ); // purposely not implemented
      /** The private copy constructor. */
      void operator=( const Self & );              // purposely not implemented

    };

    } // end namespace elastix

#ifndef ITK_MANUAL_INSTANTIATION
#include "elxTransformDirectionalStrainEnergyPenalty.hxx"
#endif

#endif // end #ifndef __elxTransformDirectionalStrainEnergyPenaltyTerm_H__


    template<class TElastix>
    class TransformDirectionalStrainEnergyPenalty :
            public itk::TransformDirectionalStrainEnergyPenaltyTerm<
                    typename MetricBase<TElastix>::FixedImageType,
                    double>,
            public MetricBase<TElastix> {
    public:

        /** Standard ITK-stuff. */
        typedef TransformDirectionalStrainEnergyPenalty Self;
        typedef itk::TransformDirectionalStrainEnergyPenaltyTerm<
                typename MetricBase<TElastix>::FixedImageType,
                double> Superclass1;
        typedef MetricBase <TElastix> Superclass2;
        typedef itk::SmartPointer <Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro( Self );

        /** Run-time type information (and related methods). */
        itkTypeMacro( TransformDirectionalStrainEnergyPenalty, itk::TransformDirectionalStrainEnergyPenaltyTerm
        );

        /** Name of this class.
         * Use this name in the parameter file to select this specific metric. \n
         * example: <tt>(Metric "TransformDirectionalStrainEnergyPenalty")</tt>\n
         */
        elxClassNameMacro( "TransformDirectionalStrainEnergyPenalty" );

        /** Typedefs from the superclass. */
        typedef typename
        Superclass1::CoordinateRepresentationType CoordinateRepresentationType;
        typedef typename Superclass1::MovingImageType MovingImageType;
        typedef typename Superclass1::MovingImagePixelType MovingImagePixelType;
        typedef typename Superclass1::MovingImageConstPointer MovingImageConstPointer;
        typedef typename Superclass1::FixedImageType FixedImageType;
        typedef typename Superclass1::FixedImageConstPointer FixedImageConstPointer;
        typedef typename Superclass1::FixedImageRegionType FixedImageRegionType;
        typedef typename Superclass1::TransformType TransformType;
        typedef typename Superclass1::TransformPointer TransformPointer;
        typedef typename Superclass1::InputPointType InputPointType;
        typedef typename Superclass1::OutputPointType OutputPointType;
        typedef typename Superclass1::TransformParametersType TransformParametersType;
        typedef typename Superclass1::TransformJacobianType TransformJacobianType;
        typedef typename Superclass1::InterpolatorType InterpolatorType;
        typedef typename Superclass1::InterpolatorPointer InterpolatorPointer;
        typedef typename Superclass1::RealType RealType;
        typedef typename Superclass1::GradientPixelType GradientPixelType;
        typedef typename Superclass1::GradientImageType GradientImageType;
        typedef typename Superclass1::GradientImagePointer GradientImagePointer;
        typedef typename Superclass1::GradientImageFilterType GradientImageFilterType;
        typedef typename Superclass1::GradientImageFilterPointer GradientImageFilterPointer;
        typedef typename Superclass1::FixedImageMaskType FixedImageMaskType;
        typedef typename Superclass1::FixedImageMaskPointer FixedImageMaskPointer;
        typedef typename Superclass1::MovingImageMaskType MovingImageMaskType;
        typedef typename Superclass1::MovingImageMaskPointer MovingImageMaskPointer;
        typedef typename Superclass1::MeasureType MeasureType;
        typedef typename Superclass1::DerivativeType DerivativeType;
        typedef typename Superclass1::ParametersType ParametersType;
        typedef typename Superclass1::FixedImagePixelType FixedImagePixelType;
        typedef typename Superclass1::MovingImageRegionType MovingImageRegionType;
        typedef typename Superclass1::ImageSamplerType ImageSamplerType;
        typedef typename Superclass1::ImageSamplerPointer ImageSamplerPointer;
        typedef typename Superclass1::ImageSampleContainerType ImageSampleContainerType;
        typedef typename
        Superclass1::ImageSampleContainerPointer ImageSampleContainerPointer;
        typedef typename Superclass1::FixedImageLimiterType FixedImageLimiterType;
        typedef typename Superclass1::MovingImageLimiterType MovingImageLimiterType;
        typedef typename
        Superclass1::FixedImageLimiterOutputType FixedImageLimiterOutputType;
        typedef typename
        Superclass1::MovingImageLimiterOutputType MovingImageLimiterOutputType;

        /** The fixed image dimension. */
        itkStaticConstMacro( FixedImageDimension,
        unsigned int,
        FixedImageType::ImageDimension );

        /** The moving image dimension. */
        itkStaticConstMacro( MovingImageDimension,
        unsigned int,
        MovingImageType::ImageDimension );

        /** Typedef's inherited from elastix. */
        typedef typename Superclass2::ElastixType ElastixType;
        typedef typename Superclass2::ElastixPointer ElastixPointer;
        typedef typename Superclass2::ConfigurationType ConfigurationType;
        typedef typename Superclass2::ConfigurationPointer ConfigurationPointer;
        typedef typename Superclass2::RegistrationType RegistrationType;
        typedef typename Superclass2::RegistrationPointer RegistrationPointer;
        typedef typename Superclass2::ITKBaseType ITKBaseType;

        /** Sets up a timer to measure the initialization time and
         * calls the Superclass' implementation.
         */
        virtual void Initialize(void) throw(itk::ExceptionObject);

        /**
         * Do some things before each resolution:
         * \li Set options for SelfHessian
         */
        virtual void BeforeEachResolution(void);

        virtual void BeforeRegistration(void);

        /**
         * Do some things after each iteration:
         * \li Print the OC, PC, LC parts of the rigidity term.
         */
        virtual void AfterEachIteration(void);

        /** This metric is advanced (so it has a sampling possibility), but it
         * purposely does not use samplers. The MetricBase class, however, issues
         * a warning if this is the case, so we overwrite that function.
         */
        virtual void SelectNewSamples(void) {}

    protected:

        /** The constructor. */
        TransformDirectionalStrainEnergyPenalty() {}

        /** The destructor. */
        virtual ~TransformDirectionalStrainEnergyPenalty() {}

    private:

        /** The private constructor. */
        TransformDirectionalStrainEnergyPenalty(const Self &); // purposely not implemented
        /** The private copy constructor. */
        void operator=(const Self &);              // purposely not implemented

    };

} // end namespace elastix

#ifndef ITK_MANUAL_INSTANTIATION

#include "elxTransformDirectionalStrainEnergyPenalty.hxx"

#endif

#endif // end #ifndef __elxTransformDirectionalStrainEnergyPenaltyTerm_H__

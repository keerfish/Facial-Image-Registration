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

#ifndef __elxTransformDirectionalStrainEnergyPenaltyTerm_HXX__
#define __elxTransformDirectionalStrainEnergyPenaltyTerm_HXX__

#include <DirectionalStrainEnergyPenalty/elxTransformDirectionalStrainEnergyPenalty.h>
#include "itkChangeInformationImageFilter.h"
#include "itkTimeProbe.h"


namespace elastix {

/**
 * ******************* BeforeRegistration ***********************
 */

    template<class TElastix>
    void
    TransformDirectionalStrainEnergyPenalty<TElastix>
    ::BeforeRegistration(void) {

        /** Read the fixed rigidity image if desired. */

        std::string fixedRigidityImageName = "";
        this->GetConfiguration()->ReadParameter(fixedRigidityImageName,
                                                "FixedRigidityImageName", this->GetComponentLabel(), 0, -1, false);

        typedef typename Superclass1::RigidityImageType RigidityImageType;
        typedef itk::ImageFileReader <RigidityImageType> RigidityImageReaderType;
        typename RigidityImageReaderType::Pointer fixedRigidityReader;
        typedef itk::ChangeInformationImageFilter <RigidityImageType> ChangeInfoFilterType;
        typedef typename ChangeInfoFilterType::Pointer ChangeInfoFilterPointer;
        typedef typename RigidityImageType::DirectionType DirectionType;

        if (fixedRigidityImageName != "") {
            /** Use the FixedRigidityImage. */
            this->SetUseFixedRigidityImage(true);

            /** Create the reader and set the filename. */
            fixedRigidityReader = RigidityImageReaderType::New();
            fixedRigidityReader->SetFileName(fixedRigidityImageName.c_str());

            /** Possibly overrule the direction cosines. */
            ChangeInfoFilterPointer infoChanger = ChangeInfoFilterType::New();
            DirectionType direction;
            direction.SetIdentity();
            infoChanger->SetOutputDirection(direction);
            infoChanger->SetChangeDirection(!this->GetElastix()->GetUseDirectionCosines());
            infoChanger->SetInput(fixedRigidityReader->GetOutput());

            /** Do the reading. */
            try {
                infoChanger->Update();
            }
            catch (itk::ExceptionObject &excp) {
                /** Add information to the exception. */
                excp.SetLocation("MattesMutualInformationWithDirectionalStrainEnergyPenalty - BeforeRegistration()");
                std::string err_str = excp.GetDescription();
                err_str += "\nError occurred while reading the fixed rigidity image.\n";
                excp.SetDescription(err_str);
                /** Pass the exception to an higher level. */
                throw excp;
            }

            /** Set the fixed rigidity image into the superclass. */
            this->SetFixedRigidityImage(infoChanger->GetOutput());
        } else {
            this->SetUseFixedRigidityImage(false);
        }

        /** Read the moving rigidity image if desired. */

        std::string movingRigidityImageName = "";
        this->GetConfiguration()->ReadParameter(movingRigidityImageName,
                                                "MovingRigidityImageName", this->GetComponentLabel(), 0, -1, false);

        typename RigidityImageReaderType::Pointer movingRigidityReader;
        if (movingRigidityImageName != "") {
            /** Use the movingRigidityImage. */
            this->SetUseMovingRigidityImage(true);

            /** Create the reader and set the filename. */
            movingRigidityReader = RigidityImageReaderType::New();
            movingRigidityReader->SetFileName(movingRigidityImageName.c_str());

            /** Possibly overrule the direction cosines. */
            ChangeInfoFilterPointer infoChanger = ChangeInfoFilterType::New();
            DirectionType direction;
            direction.SetIdentity();
            infoChanger->SetOutputDirection(direction);
            infoChanger->SetChangeDirection(!this->GetElastix()->GetUseDirectionCosines());
            infoChanger->SetInput(movingRigidityReader->GetOutput());

            /** Do the reading. */
            try {
                infoChanger->Update();
            }
            catch (itk::ExceptionObject &excp) {
                /** Add information to the exception. */
                excp.SetLocation("MattesMutualInformationWithDirectionalStrainEnergyPenalty - BeforeRegistration()");
                std::string err_str = excp.GetDescription();
                err_str += "\nError occurred while reading the moving rigidity image.\n";
                excp.SetDescription(err_str);
                /** Pass the exception to an higher level. */
                throw excp;
            }

            /** Set the moving rigidity image into the superclass. */
            this->SetMovingRigidityImage(infoChanger->GetOutput());
        } else {
            this->SetUseMovingRigidityImage(false);
        }

        /** Read the horizontalDirectionVector if desired. */

        std::string horizontalDirectionVectorImageName = "";
        this->GetConfiguration()->ReadParameter(horizontalDirectionVectorImageName,
                                                "HorizontalDirectionVectorImageName", this->GetComponentLabel(), 0, -1,
                                                false);

        typedef typename Superclass1::VectorImageType1 VectorImageType1;

        typedef itk::ImageFileReader <VectorImageType1> VectorImageReaderType1;
        typename VectorImageReaderType1::Pointer horizontalDirectionVectorReader;

        typedef itk::ChangeInformationImageFilter <VectorImageType1> ChangeInfoFilterType1;
        typedef typename ChangeInfoFilterType1::Pointer ChangeInfoFilterPointer1;
        typedef typename VectorImageType1::DirectionType DirectionType1;


        if (horizontalDirectionVectorImageName != "") {

            /** Use the horizontalDirectionVector. */
            this->SetUseHorizontalDirectionVectorImage(true);

            /** Create the reader and set the filename. */
            horizontalDirectionVectorReader = VectorImageReaderType1::New();
            horizontalDirectionVectorReader->SetFileName(horizontalDirectionVectorImageName.c_str());

            /** Possibly overrule the direction cosines. */
            ChangeInfoFilterPointer1 infoChanger = ChangeInfoFilterType1::New();
            DirectionType1 direction;
            direction.SetIdentity();
            infoChanger->SetOutputDirection(direction);
            infoChanger->SetChangeDirection(!this->GetElastix()->GetUseDirectionCosines());
            infoChanger->SetInput(horizontalDirectionVectorReader->GetOutput());

            /** Do the reading. */
            try {
                infoChanger->Update();
            }
            catch (itk::ExceptionObject &excp) {
                /** Add information to the exception. */
                excp.SetLocation("MattesMutualInformationWithDirectionalStrainEnergyPenalty - BeforeRegistration()");
                std::string err_str = excp.GetDescription();
                err_str += "\nError occurred while reading the horizontal direction vector image.\n";
                excp.SetDescription(err_str);
                /** Pass the exception to an higher level. */
                throw excp;
            }

            /** Set the horizontalDirectionVector into the superclass. */
            this->SetHorizontalDirectionVectorImage(infoChanger->GetOutput());
        } else {
            this->SetUseHorizontalDirectionVectorImage(false);
        }

        /** Read the verticalDirectionVector image if desired. */

        std::string verticalDirectionVectorImageName = "";
        this->GetConfiguration()->ReadParameter(verticalDirectionVectorImageName,
                                                "VerticalDirectionVectorImageName", this->GetComponentLabel(), 0, -1,
                                                false);


        typedef typename Superclass1::VectorImageType2 VectorImageType2;
        typedef itk::ImageFileReader <VectorImageType2> VectorImageReaderType2;
        typename VectorImageReaderType2::Pointer verticalDirectionVectorReader;

        typedef itk::ChangeInformationImageFilter <VectorImageType2> ChangeInfoFilterType2;
        typedef typename ChangeInfoFilterType2::Pointer ChangeInfoFilterPointer2;
        typedef typename VectorImageType2::DirectionType DirectionType2;


        if (verticalDirectionVectorImageName != "") {
            /** Use the VerticalDirectionVectorImage. */
            this->SetUseVerticalDirectionVectorImage(true);

            /** Create the reader and set the filename. */
            verticalDirectionVectorReader = VectorImageReaderType2::New();
            verticalDirectionVectorReader->SetFileName(verticalDirectionVectorImageName.c_str());

            /** Possibly overrule the direction cosines. */
            ChangeInfoFilterPointer2 infoChanger = ChangeInfoFilterType2::New();
            DirectionType2 direction;
            direction.SetIdentity();
            infoChanger->SetOutputDirection(direction);
            infoChanger->SetChangeDirection(!this->GetElastix()->GetUseDirectionCosines());
            infoChanger->SetInput(verticalDirectionVectorReader->GetOutput());

            /** Do the reading. */
            try {
                infoChanger->Update();
            }
            catch (itk::ExceptionObject &excp) {
                /** Add information to the exception. */
                excp.SetLocation("MattesMutualInformationWithDirectionalStrainEnergyPenalty - BeforeRegistration()");
                std::string err_str = excp.GetDescription();
                err_str += "\nError occurred while reading the vertical direction vector image.\n";
                excp.SetDescription(err_str);
                /** Pass the exception to an higher level. */
                throw excp;
            }

            /** Set the VerticalDirectionVectorImage into the superclass. */
            this->SetVerticalDirectionVectorImage(infoChanger->GetOutput());
        } else {
            this->SetUseVerticalDirectionVectorImage(false);
        }



        /** Important check: at least one image must be given. */
        if (fixedRigidityImageName == "" && movingRigidityImageName == "" && horizontalDirectionVectorImageName == "" &&
            verticalDirectionVectorImageName == "") {
            xl::xout["warning"] << "WARNING: FixedRigidityImageName and "
                                << "MovingRigidityImage are both not supplied.\n"
                                << "VerticalDirectionVectorImage are both not supplied.\n"
                                << "  The strain energy penalty term is evaluated on entire input "
                                << "transform domain." << std::endl;
        }

        /** Add target cells to xout["iteration"]. */

    } // end BeforeRegistration()


/**
 * ******************* Initialize ***********************
 */

    template<class TElastix>
    void
    TransformDirectionalStrainEnergyPenalty<TElastix>
    ::Initialize(void) throw(itk::ExceptionObject) {
        itk::TimeProbe timer;
        timer.Start();
        this->Superclass1::Initialize();
        timer.Stop();
        elxout << "Initialization of TransformDirectionalStrainEnergyPenalty metric took: "
               << static_cast< long >( timer.GetMean() * 1000 ) << " ms." << std::endl;

        /** Check stuff. */
        this->CheckUseAndCalculationBooleans();

    } // end Initialize()


/**
 * ***************** BeforeEachResolution ***********************
 */

    template<class TElastix>
    void
    TransformDirectionalStrainEnergyPenalty<TElastix>
    ::BeforeEachResolution(void) {
        /** Get the current resolution level. */
        unsigned int level
                = this->m_Registration->GetAsITKBaseType()->GetCurrentLevel();

        /** Get and set the dilateRigidityImages. */
        bool dilateRigidityImages = true;
        this->GetConfiguration()->ReadParameter(dilateRigidityImages,
                                                "DilateRigidityImages", this->GetComponentLabel(), level, 0);
        this->SetDilateRigidityImages(dilateRigidityImages);

        bool dilateVectorImages = true;
        this->GetConfiguration()->ReadParameter(dilateVectorImages,
                                                "DilateVectorImages", this->GetComponentLabel(), level, 0);
        this->SetDilateVectorImages(dilateVectorImages);

        /** Get and set the dilationRadiusMultiplier. */
        double dilationRadiusMultiplier = 1.0;
        this->GetConfiguration()->ReadParameter(dilationRadiusMultiplier,
                                                "DilationRadiusMultiplier", this->GetComponentLabel(), level, 0);
        this->SetDilationRadiusMultiplier(dilationRadiusMultiplier);

        /** Get and set the usage of the GeneralStrain condition part. */
        bool useGeneralStrainPart1 = true;
        this->GetConfiguration()->ReadParameter(useGeneralStrainPart1,
                                                "UseGeneralStrainPart1", this->GetComponentLabel(), level, 0);
        this->SetUseGeneralStrainPart1(useGeneralStrainPart1);

        bool useGeneralStrainPart2 = true;
        this->GetConfiguration()->ReadParameter(useGeneralStrainPart2,
                                                "UseGeneralStrainPart2", this->GetComponentLabel(), level, 0);
        this->SetUseGeneralStrainPart2(useGeneralStrainPart2);

        /** Set the usage of the PrincipleStrain condition part. */
        bool usePrincipleStrainPart1 = true;
        this->GetConfiguration()->ReadParameter(usePrincipleStrainPart1,
                                                "UsePrincipleStrainPart1", this->GetComponentLabel(), level, 0);
        this->SetUsePrincipleStrainPart1(usePrincipleStrainPart1);

        bool usePrincipleStrainPart2 = true;
        this->GetConfiguration()->ReadParameter(usePrincipleStrainPart2,
                                                "UsePrincipleStrainPart2", this->GetComponentLabel(), level, 0);
        this->SetUsePrincipleStrainPart2(usePrincipleStrainPart2);

        bool usePrincipleStrainPart3 = true;
        this->GetConfiguration()->ReadParameter(usePrincipleStrainPart3,
                                                "UsePrincipleStrainPart3", this->GetComponentLabel(), level, 0);
        this->SetUsePrincipleStrainPart3(usePrincipleStrainPart3);

        /** Set the calculation of the GeneralStrain condition part. */
        bool calculateGeneralStrainPart1 = true;
        this->GetConfiguration()->ReadParameter(calculateGeneralStrainPart1,
                                                "CalculateGeneralStrainPart1", this->GetComponentLabel(), level, 0);
        this->SetCalculateGeneralStrainPart1(calculateGeneralStrainPart1);

        bool calculateGeneralStrainPart2 = true;
        this->GetConfiguration()->ReadParameter(calculateGeneralStrainPart2,
                                                "CalculateGeneralStrainPart2", this->GetComponentLabel(), level, 0);
        this->SetCalculateGeneralStrainPart2(calculateGeneralStrainPart2);

        /** Set the calculation of the PrincipleStrain condition part. */
        bool calculatePrincipleStrainPart1 = true;
        this->GetConfiguration()->ReadParameter(calculatePrincipleStrainPart1,
                                                "CalculatePrincipleStrainPart1", this->GetComponentLabel(), level, 0);
        this->SetCalculatePrincipleStrainPart1(calculatePrincipleStrainPart1);

        bool calculatePrincipleStrainPart2 = true;
        this->GetConfiguration()->ReadParameter(calculatePrincipleStrainPart2,
                                                "CalculatePrincipleStrainPart2", this->GetComponentLabel(), level, 0);
        this->SetCalculatePrincipleStrainPart2(calculatePrincipleStrainPart2);

        bool calculatePrincipleStrainPart3 = true;
        this->GetConfiguration()->ReadParameter(calculatePrincipleStrainPart3,
                                                "CalculatePrincipleStrainPart3", this->GetComponentLabel(), level, 0);
        this->SetCalculatePrincipleStrainPart3(calculatePrincipleStrainPart3);


        double preSetParameter = 0.0;
        this->m_Configuration->ReadParameter(preSetParameter,
                                             "PreSetParameter", this->GetComponentLabel(), level, 0);
        this->SetPreSetParameter(preSetParameter);

        double generalStrainPartWeight = 1.0;
        this->m_Configuration->ReadParameter(generalStrainPartWeight,
                                             "GeneralStrainPartWeight", this->GetComponentLabel(), level, 0);
        this->SetGeneralStrainPartWeight(generalStrainPartWeight);

        /** Set the PrincipleStrainWeight of this level. */
        double principleStrainPartWeight = 1.0;
        this->m_Configuration->ReadParameter(principleStrainPartWeight,
                                             "PrincipleStrainPartWeight", this->GetComponentLabel(), level, 0);
        this->SetPrincipleStrainPartWeight(principleStrainPartWeight);


        double horizontalDirectionWeight = 1.0;
        this->m_Configuration->ReadParameter(horizontalDirectionWeight,
                                             "horizontalDirectionWeight", this->GetComponentLabel(), level, 0);
        this->SethorizontalDirectionWeight(horizontalDirectionWeight);


        double verticalDirectionWeight = 1.0;
        this->m_Configuration->ReadParameter(verticalDirectionWeight,
                                             "verticalDirectionWeight", this->GetComponentLabel(), level, 0);
        this->SetverticalDirectionWeight(verticalDirectionWeight);

    } // end BeforeEachResolution()


/**
 * ***************AfterEachIteration ****************************
 */

    template<class TElastix>
    void
    TransformDirectionalStrainEnergyPenalty<TElastix>
    ::AfterEachIteration(void) {

    } // end AfterEachIteration()


} // end namespace elastix

#endif // end #ifndef __elxTransformDirectionalStrainEnergyPenaltyTerm_HXX__


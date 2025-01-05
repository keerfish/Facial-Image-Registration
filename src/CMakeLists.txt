cmake_minimum_required(VERSION 3.9)
project(CodeDirectionalStrainEnergy_PenaltyTerm)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        Example/fixed.png
        Example/moving.png
        Example/Parameters.txt
        elxTransformDirectionalStrainEnergyPenalty.cxx
        elxTransformDirectionalStrainEnergyPenalty.h
        elxTransformDirectionalStrainEnergyPenalty.hxx
        itkTransformDirectionalStrainEnergyPenaltyTerm.h
        itkTransformDirectionalStrainEnergyPenaltyTerm.hxx)

add_executable(CodeDirectionalStrainEnergy_PenaltyTerm ${SOURCE_FILES})
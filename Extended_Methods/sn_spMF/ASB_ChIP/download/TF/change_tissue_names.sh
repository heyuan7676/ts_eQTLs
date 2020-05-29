#!/bin/bash
## Hand curation

file="$1"

sed -i "/gastrocnemius medialis/d" ${file}
sed -i "/Peyer's patch/d" ${file}
sed -i "/occipital lobe/d" ${file}
sed -i "/globus pallidus/d" ${file}
sed -i "/pons/d" ${file}
sed -i "/Ammon's horn/d" ${file}
sed -i "/thymus/d" ${file}
sed -i "/right atrium auricular region/d" ${file}
sed -i "/arm bone/d" ${file}
sed -i "/femur/d" ${file}
sed -i "/large intestine/d" ${file}
sed -i "/intestine/d" ${file}
sed -i "/leg bone/d" ${file}
sed -i "/limb/d" ${file}
sed -i "/medulla oblongata/d" ${file}
sed -i "/superior temporal gyrus/d" ${file}
sed -i "/renal cortex interstitium/d" ${file}
sed -i "/middle frontal gyrus/d" ${file}
sed -i "/renal pelvis/d" ${file}
sed -i "/placenta/d" ${file}
sed -i "/posterior cingulate/d" ${file}
sed -i "/breast epithelium/d" ${file}
sed -i "/chorion/d" ${file}
sed -i "/retina/d" ${file}
sed -i "/umbilical cord/d" ${file}
sed -i "/germinal center/d" ${file}
sed -i "/esophagus squamous epithelium/d" ${file}
sed -i "/eye/d" ${file}
sed -i "/inferior parietal cortex/d" ${file}
sed -i "/islet of Langerhans/d" ${file}
sed -i "/tongue/d" ${file}
sed -i "/embryonic facial prominence/d" ${file}
sed -i "/spinal cord/d" ${file}
sed -i "/midBrain/d" ${file}
sed -i "/bone marrow/d" ${file}
sed -i "/cortical plate/d" ${file}
sed -i "/forebrain/d" ${file}
sed -i "/hindbrain/d" ${file}
sed -i "/midbrain/d" ${file}
sed -i "/olfactory bulb/d" ${file}
sed -i "/parathyroid adenoma/d" ${file}
sed -i "/suprapubic skin/d" ${file}


sed -i "s/right kidney/Kidney_Cortex/g" ${file}
sed -i "s/left kidney/Kidney_Cortex/g" ${file}
sed -i "s/kidney/Kidney_Cortex/g" ${file}

sed -i "s/gastroesophageal sphincter/Esophagus_Gastroesophageal_Junction/g" ${file}
sed -i "s/esophagus muscularis mucosa/Esophagus_Mucosa/g" ${file}
sed -i "s/coronary artery/Artery_Coronary/g" ${file}
sed -i "s/adrenal gland/Adrenal_Gland/g" ${file}
sed -i "s/uterus/Uterus/g" ${file}
sed -i "s/vagina/Vagina/g" ${file}
sed -i "s/ovary/Ovary/g" ${file}
sed -i "s/prostate gland/Prostate/g" ${file}
sed -i "s/prostate/Prostate/g" ${file}
sed -i "s/stomach/Stomach/g" ${file}
sed -i "s/tibial nerve/Nerve_Tibial/g" ${file}
sed -i "s/spleen/Spleen/g" ${file}
sed -i "s/thyroid gland/Thyroid/g" ${file}
sed -i "s/ascending aorta/Artery_Aorta/g" ${file}
sed -i "s/thoracic aorta/Artery_Aorta/g" ${file}
sed -i "s/tibial artery/Artery_Tibial/g" ${file}
sed -i "s/lower leg skin/Skin_Sun_Exposed_Lower_leg/g" ${file}
sed -i "s/skin of body/Skin_Sun_Exposed_Lower_leg/g" ${file}
sed -i "s/small intestine/Small_Intestine_Terminal_Ileum/g" ${file}



sed -i "s/omental fat pad/Adipose_Subcutaneous/g" ${file}
sed -i "s/subcutaneous adipose tissue/Adipose_Subcutaneous/g" ${file}


sed -i "s/heart left ventricle/Heart_Left_Ventricle/g" ${file}
sed -i "s/heart right ventricle/Heart_Left_Ventricle/g" ${file}
sed -i "s/heart/Heart_Left_Ventricle/g" ${file}
sed -i "s/left cardiac atrium/Heart_Atrial_Appendage/g" ${file}

sed -i "s/upper lobe of left lung/Lung/g" ${file}
sed -i "s/right lung/Lung/g" ${file}
sed -i "s/left lung/Lung/g" ${file}
sed -i "s/lung/Lung/g" ${file}

sed -i "s/body of pancreas/Pancreas/g" ${file}
sed -i "s/pancreas/Pancreas/g" ${file}

sed -i "s/cerebellar cortex/Brain_Cerebellar_Hemisphere/g" ${file}
sed -i "s/cerebellum/Brain_Cerebellum/g" ${file}
sed -i "s/putamen/Brain_Putamen_basal_ganglia/g" ${file}
sed -i "s/frontal cortex/Brain_Frontal_Cortex_BA9/g" ${file}
sed -i "s/brain/Brain/g" ${file}
sed -i "s/caudate nucleus/Brain_Caudate_basal_ganglia/g" ${file}


sed -i "s/muscle of arm/Muscle_Skeletal/g" ${file}
sed -i "s/psoas muscle/Muscle_Skeletal/g" ${file}
sed -i "s/muscle of back/Muscle_Skeletal/g" ${file}
sed -i "s/forelimb muscle/Muscle_Skeletal/g" ${file}
sed -i "s/hindlimb muscle/Muscle_Skeletal/g" ${file}
sed -i "s/muscle of leg/Muscle_Skeletal/g" ${file}
sed -i "s/muscle of trunk/Muscle_Skeletal/g" ${file}


sed -i "s/right lobe of liver/Liver/g" ${file}
sed -i "s/liver/Liver/g" ${file}

sed -i "s/sigmoid colon/Colon_Sigmoid/g" ${file}
sed -i "s/transverse colon/Colon_Transverse/g" ${file}
sed -i "s/testis/Testis/g" ${file}

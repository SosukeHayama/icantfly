//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
//			�A���e�i�`��w�b�_�[�@�@antenna.h
//			�ђʌ`�A���e�i
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
//				�v�Z
//////////////////////////////////////////////////////////////////////////////////
void Calculation(void)
{
	for(int AA = 59;AA<= 59; ++AA){		//x������̋���d�p(�Œ肷��Ƃ�ver)
										//d_pole = 0.001 * RAMDA0 * AA;
	for(int CC = 0; CC<= 0; ++CC){		//f���p
	for(int DD = 125; DD<= 125; ++DD){	//F�`#1�Ԃ̒������̒��œK���p
										//F_1_pole = 0.001*DD * RAMDA0;1.169
	for(int EE = 1071; EE<= 1071; ++EE){	//C�œK���p
	for(int BB = 1023; BB<= 1023; ++BB){	//Lx+C/4�œK���p
										//Lx_pole = 0.001*BB * RAMDA0;
	//========================================================
	//			���[�J���ϐ�
	//========================================================
	//--�J�E���^
	int i,j;

	//--�v�Z�Ŏg�p����ϐ�
	double rrx , rry , rrz;				//�Z�O�����g�n�_
	double drx , dry , drz;				//�Z�O�����g�I�_
	double ssx , ssy , ssz;				//�Z�O�����g�̒P�ʃx�N�g��
	double dcl;							//�Z�O�����g����

	//--�G�������g���֌W�i���C���[�{���j
	int n_array;

	//--�e�����̕�����
	int nbun_l , nbun_l_o;			//�����`��1�ӂ̕�����
	int nbun_d;						//x���Ɛ������̂̊Ԋud���̕�����
	int nbun_lps , nbun_lps_o;		//(�����`��1�ӂ̒���/2 + �Ԋud)�̕�����
	int nbun_lms , nbun_lms_o;		//(�����`��1�ӂ̒���/2 - �Ԋud)�̕�����
	int nbun_h , nbun_h_o;			//���d���������̕�����
	int nbun_F_1 , nbun_F_1_o;		//���[�v#1�Ƌ��d���̊Ԃ̃��C��
	int nbun_Lx , nbun_Lx_o;		//��������Lx
	int nbun_F;						//���d���������c


	//--�e���̃Z�O�����g��
	int n_pole;

	//--�݌v�l�֌W
	double wire_radius;			//���C���[���a
	double C_pole;				//�����`���[�v���͒�
	double l_pole;				//�����`��1�ӂ̒���
	double d_pole;				//x���Ɛ������̂̊Ԋud
	double lpd_pole;			//�����`��1�ӂ̒���/2 + �Ԋud
	double lmd_pole;			//�����`��1�ӂ̒���/2 - �Ԋud
	double h_pole;				//����h
	double F_1_pole;			//���[�v#1�Ƌ��d���̊Ԃ̃��C��
	double Lx_pole;				//��������Lx
	double F_pole;				//���d���������c

	int hensuu;					//�Ȃ�ɂł��g���Ă����ϐ�
	int hensuu2;

	int soshi_no = 4;

	//========================================================
	//			���ˊE�v�Z�ݒ�		�P�ʁF(deg)
	//========================================================
	//--�v�Z���[�h
	AXMODE = 0;			//(0=�ӌŒ�, 1=�ƌŒ�)

	//--�ώ�
	DegDelta =    1.0;	//���ݕ�("0.0"�֎~)
	DegStart =  -90.0;	//�����p
	DegWidth =  180.0;	//�͈�

	//--�Œ莲
	FixAngle =  0.0;	//�Œ�p�x

	//========================================================
	//			�݌v���g��		FREQ0	(f0)	�P�ʁF(Hz)
	//			���d���g��		USEF	(f)		�P�ʁF(Hz)
	//========================================================
	//--�݌v���g���ݒ�@
	FREQ0 = 3.081 * pow(10.0,9.0);

	//--�d�������ݒ�
	//USEF  = FREQ0 * 1.0;					//f���ȊO�̂Ƃ�
	USEF  = FREQ0 * (1.000+0.001*CC);					//f���̂Ƃ�

	//--�݌v���g���̎��R��Ԕg��
	RAMDA0= C/FREQ0;

	//========================================================
	//			�`��l��`
	//			������RAMDA0����@�W���Œl��`(������)
	//========================================================

	C_pole = 0.001 * EE * RAMDA0;				//�����`���[�v���͒�
	l_pole = C_pole/4.0;				//�����`��1�ӂ̒���

	d_pole = 0.001 * RAMDA0 * AA;			//x���Ɛ������̂̊Ԋud
	lpd_pole = l_pole/2.0 + d_pole;		//�����`��1�ӂ̒���/2 + �Ԋud
	lmd_pole = l_pole/2.0 - d_pole;		//�����`��1�ӂ̒���/2 - �Ԋud

	h_pole = 0.125 * RAMDA0;			//����h

	F_1_pole = 0.001*DD * RAMDA0;			//���[�v#1�Ƌ��d���̊Ԃ̃��C��
	Lx_pole = 0.001*BB * RAMDA0 - C_pole/4.0;			//��������Lx
	F_pole = 0.125 * RAMDA0;			//���d���������c

	wire_radius = 0.005 * RAMDA0;	//���C���[���a

	//========================================================
	//			�o�̓t�@�C���p�̃p�����[�^
	//			PARA1 �� PARA2 �ɒl����͂���D
	//========================================================
	PARA1 = d_pole/RAMDA0;
	PARA2 = C_pole / RAMDA0;
	PARA3 = (Lx_pole+l_pole)/RAMDA0;
	
	//PARA1 = CC;
	PARA4 = USEF/FREQ0;


	//========================================================
	//			�S���C���[��	NWIR
	//			�Z�O�����g��	NSEG
	//			���d�_��		NFED
	//			�S�d���v�Z�_��	NSEG0
	//			�y���א�		NLOAD		
	//========================================================

	//--�S���C���[���ݒ�
	n_array = (2*soshi_no					//���[�v
			+ 1*soshi_no					//��������
			+ 1*(soshi_no-1))*2				//��������Lx
			+ 1;							//���d��
	NWIR = n_array;
  	double segtyou=0.025*RAMDA0;
	
	//--����������
	nbun_l = int(l_pole/segtyou+0.5);							//�����`��1�ӂ̕�����
	nbun_d = int(d_pole/segtyou+0.5);								//x���Ɛ������̂̊Ԋud���̕�����
	nbun_lps = int(lpd_pole/segtyou+0.5);			//(�����`��1�ӂ̒���/2 + �Ԋud)�̕�����
	nbun_lms = int(lmd_pole/segtyou+0.5);			//(�����`��1�ӂ̒���/2 - �Ԋud)�̕�����
	nbun_h = int(h_pole/segtyou+0.5);								//��������F-F'�̕�����
	nbun_F_1 = int(F_1_pole/segtyou+0.5);							//���[�v#1�Ƌ��d���̊Ԃ̃��C��
	nbun_Lx = int(Lx_pole/segtyou+0.5);							//��������Lx
	nbun_F = int(F_pole/segtyou+0.5);								//���d���������c

	//--�I�[�o�[���b�v��������
	nbun_l_o = nbun_l + 2;					//�����`��1�ӂ̕�����
	nbun_lps_o = nbun_lps + 1;				//(�����`��1�ӂ̒���/2 + �Ԋud)�̕�����
	nbun_lms_o = nbun_lms + 1;				//(�����`��1�ӂ̒���/2 - �Ԋud)�̕�����
	nbun_h_o = nbun_h + 1;					//��������F-F'�̕�����
	nbun_F_1_o = nbun_F_1 + 1;				//���[�v#1�Ƌ��d���̊Ԃ̃��C��
	nbun_Lx_o = nbun_Lx + 2;				//��������Lx


	//--�S�Z�O�����g���ݒ�
	n_pole = ( (nbun_lps_o + nbun_l + nbun_lps + nbun_lms_o + nbun_l + nbun_lms)*soshi_no	//���[�v
			+ nbun_l_o*soshi_no																//��������
			+ nbun_Lx_o*(soshi_no-1) )*2													//��������Lx
			+ (nbun_F_1_o + nbun_h)*2;

	NSEG = n_pole;						//�S�Z�O�����g��

	//--�S�d���v�Z�_��  [NSEG0=NSEG-NWIR]	�������������폜�@�s��
	NSEG0 = NSEG - NWIR;

	//--���d�_��
	NFED = 1;

	//--�C���s�[�_���X���א�
	NLOAD = 0;

	//========================================================
	//			�z��m�ۂƏ�����
	//			[���ˊE�v�Z�ݒ�,NWIR,NSEG,NFED,NLOAD]
	//				�ݒ��Ɋ֐��Ăяo��
	//			�������������폜�@�s��
	//========================================================
	MakeMatAll();								//�z��m��(���R�����g�A�E�g�֎~��)
	Initialization();							//�z��ƌv�Z�p�ϐ��̏�����
	
	//========================================================
	//			�݌v�l����		�P�ʁF(m)
	//			�Z�O�����g�n�_	RX[ ] , RY[ ] , RZ[ ] 
	//			�P�ʃx�N�g��	SX[ ] , SY[ ] , SZ[ ]
	//			�Z�O�����g��	SEGL[ ]
	//========================================================

	////���˔\////
	//���˔�_�\_���[�v_�㕔_�I�[�o�[���b�v
		i = 0;
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = F_1_pole;							//�n�_�ʒu�v�Zx
		rry = -(lmd_pole/nbun_lms) - d_pole;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = F_1_pole;							//�n�_�ʒu�v�Zx
		dry = -d_pole;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��

	//���˔�_�\_���[�v_�㕔_����
	for(j = 0; j < nbun_lps_o - 1; ++j){
		i = j + (1);
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = F_1_pole;							//�n�_�ʒu�v�Zx
		rry = (lpd_pole/nbun_lps)*(j+0) - d_pole;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = F_1_pole;							//�n�_�ʒu�v�Zx
		dry = (lpd_pole/nbun_lps)*(j+1) - d_pole;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
	}

	//���˔�_�\_���[�v_�㕔_���
	for(j = 0; j < nbun_l; ++j){
		i = j + (nbun_lps_o);
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = (l_pole/nbun_l)*(j+0) + F_1_pole;							//�n�_�ʒu�v�Zx
		rry = l_pole/2.0;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = (l_pole/nbun_l)*(j+1) + F_1_pole;							//�n�_�ʒu�v�Zx
		dry = l_pole/2.0;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
	}

	//���˔�_�\_���[�v_�㕔_�E��
	for(j = 0; j < nbun_lps; ++j){
		i = j + (nbun_lps_o + nbun_l);
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = F_1_pole + l_pole;							//�n�_�ʒu�v�Zx
		rry = -(lpd_pole/nbun_lps)*(j+0) + l_pole/2.0;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = F_1_pole + l_pole;							//�n�_�ʒu�v�Zx
		dry = -(lpd_pole/nbun_lps)*(j+1) + l_pole/2.0;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
	}

	//���˔�_�\_���[�v_����_�I�[�o�[���b�v
		i = (nbun_lps_o + nbun_l + nbun_lps);
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = F_1_pole + l_pole;							//�n�_�ʒu�v�Zx
		rry = (lpd_pole/nbun_lps) - d_pole;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = F_1_pole + l_pole;							//�n�_�ʒu�v�Zx
		dry = -d_pole;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��

	//���˔�_�\_���[�v_����_�E��
	for(j = 0; j < nbun_lms_o - 1; ++j){
		i = j + (nbun_lps_o + nbun_l + nbun_lps) + (1);
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = F_1_pole + l_pole;							//�n�_�ʒu�v�Zx
		rry = - (lmd_pole/nbun_lms)*(j+0) -d_pole;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = F_1_pole + l_pole;							//�n�_�ʒu�v�Zx
		dry = - (lmd_pole/nbun_lms)*(j+1) -d_pole;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
	}
	
	//���˔�_�\_���[�v_����_����
	for(j = 0; j < nbun_l; ++j){
		i = j + (nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o);
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = -(l_pole/nbun_l)*(j+0) + F_1_pole + l_pole;							//�n�_�ʒu�v�Zx
		rry = -l_pole/2.0;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = -(l_pole/nbun_l)*(j+1) + F_1_pole + l_pole;							//�n�_�ʒu�v�Zx
		dry = -l_pole/2.0;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
	}

	//���˔�_�\_���[�v_����_����
	for(j = 0; j < nbun_lms; ++j){
		i = j + (nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l);
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = F_1_pole;							//�n�_�ʒu�v�Zx
		rry = (lmd_pole/nbun_lms)*(j+0) -l_pole/2.0;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = F_1_pole;							//�n�_�ʒu�v�Zx
		dry = (lmd_pole/nbun_lms)*(j+1) -l_pole/2.0;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
	}

	//���˔�_�\_���[�v#N
	for(hensuu = 2; hensuu <= soshi_no; ++hensuu){
		for(j = 0; j < (nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms); ++j){
			i = j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*(hensuu-1);
			//--�v�Z���ʑ��
			RX[i] = RX[j] + (l_pole + Lx_pole)*(hensuu-1);				//�n�_�ʒux
			RY[i] = RY[j];											//�n�_�ʒuy
			RZ[i] = RZ[j];											//�n�_�ʒuz
			SX[i] = SX[j];											//�P�ʃx�N�g��x
			SY[i] = SY[j];											//�P�ʃx�N�g��y
			SZ[i] = SZ[j];											//�P�ʃx�N�g��z
			SEGL[i] = SEGL[j];										//�Z�O�����g��
		}
	}

	//���˔�_�\_��������_�I�[�o�[���b�v
		i = ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no;
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = F_1_pole;							//�n�_�ʒu�v�Zx
		rry = (lpd_pole/nbun_lps)-d_pole;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = F_1_pole;							//�n�_�ʒu�v�Zx
		dry = -d_pole;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��

	//���˔�_�\_��������
	for(j = 0; j < nbun_l; ++j){
		i = j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ 1;
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = (l_pole/nbun_l)*(j+0) + F_1_pole;							//�n�_�ʒu�v�Zx
		rry = -d_pole;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx =(l_pole/nbun_l)*(j+1) + F_1_pole;								//�n�_�ʒu�v�Zx
		dry = -d_pole;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
	}

	//���˔�_�\_��������_�I�[�o�[���b�v
		i = ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l + 1;
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = F_1_pole + l_pole;							//�n�_�ʒu�v�Zx
		rry = -d_pole;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = F_1_pole + l_pole;							//�n�_�ʒu�v�Zx
		dry = -(lmd_pole/nbun_lms)-d_pole;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��

	//���˔�_�\_��������_#2�`#N
	for(hensuu = 2; hensuu <= soshi_no; ++hensuu){
		for(j = 0; j < nbun_l_o; ++j){
			i = j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
				+ nbun_l_o*(hensuu-1);
			//--�v�Z���ʑ��
			RX[i] = RX[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no] + (l_pole + Lx_pole)*(hensuu-1);				//�n�_�ʒux
			RY[i] = RY[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no];											//�n�_�ʒuy
			RZ[i] = RZ[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no];											//�n�_�ʒuz
			SX[i] = SX[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no];											//�P�ʃx�N�g��x
			SY[i] = SY[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no];											//�P�ʃx�N�g��y
			SZ[i] = SZ[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no];											//�P�ʃx�N�g��z
			SEGL[i] = SEGL[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no];										//�Z�O�����g��
		}
	}

	//���˔�_�\_��������Lx_�I�[�o�[���b�v
		i = ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no;
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = -(l_pole/nbun_l) + F_1_pole + l_pole;							//�n�_�ʒu�v�Zx
		rry = -d_pole;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = F_1_pole + l_pole;							//�n�_�ʒu�v�Zx
		dry = -d_pole;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��

	//���˔�_�\_��������Lx
	for(j = 0; j < nbun_Lx; ++j){
		i = j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ 1;
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = (Lx_pole/nbun_Lx)*(j+0) + F_1_pole + l_pole;							//�n�_�ʒu�v�Zx
		rry = -d_pole;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx =(Lx_pole/nbun_Lx)*(j+1) + F_1_pole + l_pole;								//�n�_�ʒu�v�Zx
		dry = -d_pole;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
	}

	//���˔�_�\_��������Lx_�I�[�o�[���b�v
		i = ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx + 1;
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = F_1_pole + l_pole + Lx_pole;							//�n�_�ʒu�v�Zx
		rry = -d_pole;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = (l_pole/nbun_l) + F_1_pole + l_pole + Lx_pole;							//�n�_�ʒu�v�Zx
		dry = -d_pole;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��

	//���˔�_�\_��������Lx_�c��
	for(hensuu = 2; hensuu <= soshi_no-1; ++hensuu){
		for(j = 0; j < nbun_Lx_o; ++j){
		i = j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx_o*(hensuu-1);
			//--�v�Z���ʑ��
			RX[i] = RX[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no] + (l_pole + Lx_pole)*(hensuu-1);				//�n�_�ʒux
			RY[i] = RY[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no];											//�n�_�ʒuy
			RZ[i] = RZ[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no];											//�n�_�ʒuz
			SX[i] = SX[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no];											//�P�ʃx�N�g��x
			SY[i] = SY[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no];											//�P�ʃx�N�g��y
			SZ[i] = SZ[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no];											//�P�ʃx�N�g��z
			SEGL[i] = SEGL[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no];										//�Z�O�����g��
		}
	}

	//���˔�
		for(j = 0; j < ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
						+ nbun_l_o*soshi_no
						+ nbun_Lx_o*(soshi_no-1); ++j){
			i = j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
				+ nbun_l_o*soshi_no
				+ nbun_Lx_o*(soshi_no-1);
			//--�v�Z���ʑ��
			RX[i] = RX[j];											//�n�_�ʒux
			RY[i] = RY[j];											//�n�_�ʒuy
			RZ[i] = -RZ[j];											//�n�_�ʒuz
			SX[i] = SX[j];											//�P�ʃx�N�g��x
			SY[i] = SY[j];											//�P�ʃx�N�g��y
			SZ[i] = SZ[j];											//�P�ʃx�N�g��z
			SEGL[i] = SEGL[j];										//�Z�O�����g��
		}

	//���d��_�I�[�o�[���b�v
		i = ( ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx_o*(soshi_no-1) )*2;
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = (l_pole/nbun_l) + F_1_pole;							//�n�_�ʒu�v�Zx
		rry = -d_pole;							//�n�_�ʒu�v�Zy
		rrz = -h_pole;														//�n�_�ʒu�v�Zz
		drx = F_1_pole;							//�n�_�ʒu�v�Zx
		dry = -d_pole;							//�n�_�ʒu�v�Zy
		drz = -h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��

	//���d��_F�`#1_���˔�
	for(j = 0; j < nbun_F_1; ++j){
		i = j + ( ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx_o*(soshi_no-1) )*2
			+ 1;
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = -(F_1_pole/nbun_F_1)*(j+0) + F_1_pole;							//�n�_�ʒu�v�Zx
		rry = -d_pole;							//�n�_�ʒu�v�Zy
		rrz = -h_pole;														//�n�_�ʒu�v�Zz
		drx = -(F_1_pole/nbun_F_1)*(j+1) + F_1_pole;								//�n�_�ʒu�v�Zx
		dry = -d_pole;							//�n�_�ʒu�v�Zy
		drz = -h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
	}

	//���d��_������
	for(j = 0; j < nbun_h*2; ++j){
		i = j + ( ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx_o*(soshi_no-1) )*2
			+ nbun_F_1_o;
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = 0;							//�n�_�ʒu�v�Zx
		rry = - d_pole;							//�n�_�ʒu�v�Zy
		rrz = (h_pole/nbun_h)*(j+0)-h_pole;														//�n�_�ʒu�v�Zz
		drx = 0;								//�n�_�ʒu�v�Zx
		dry = - d_pole;							//�n�_�ʒu�v�Zy
		drz = (h_pole/nbun_h)*(j+1)-h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
	}

	//���d��_F�`#1_���˔\
	for(j = 0; j < nbun_F_1; ++j){
		i = j + ( ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx_o*(soshi_no-1) )*2
			+ nbun_F_1_o + nbun_h*2;
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = (F_1_pole/nbun_F_1)*(j+0);							//�n�_�ʒu�v�Zx
		rry = -d_pole;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = (F_1_pole/nbun_F_1)*(j+1);								//�n�_�ʒu�v�Zx
		dry = -d_pole;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
	}

	//���d��_�I�[�o�[���b�v
		i = ( ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx_o*(soshi_no-1) )*2
			+ nbun_F_1_o + nbun_h*2 + nbun_F_1;
		//--�n�_ �I�_ �Z�O�����g���v�Z
		rrx = F_1_pole;							//�n�_�ʒu�v�Zx
		rry = -d_pole;							//�n�_�ʒu�v�Zy
		rrz = h_pole;														//�n�_�ʒu�v�Zz
		drx = (l_pole/nbun_l) + F_1_pole;							//�n�_�ʒu�v�Zx
		dry = -d_pole;							//�n�_�ʒu�v�Zy
		drz = h_pole;														//�I�_�ʒu�v�Zz
		//--�Z�O�����g�� �P�ʕ����x�N�g���̊e�����̌v�Z
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//�Z�O�����g���v�Z
		ssx = (drx-rrx) / dcl;									//�P�ʃx�N�g��x
		ssy = (dry-rry) / dcl;									//�P�ʃx�N�g��y
		ssz = (drz-rrz) / dcl;									//�P�ʃx�N�g��z
		//--�v�Z���ʑ��
		RX[i] = rrx;											//�n�_�ʒux
		RY[i] = rry;											//�n�_�ʒuy
		RZ[i] = rrz;											//�n�_�ʒuz
		SX[i] = ssx;											//�P�ʃx�N�g��x
		SY[i] = ssy;											//�P�ʃx�N�g��y
		SZ[i] = ssz;											//�P�ʃx�N�g��z
		SEGL[i] = dcl;											//�Z�O�����g��
	
	//========================================================
	//			�e�Z�O�����g�̃��C���[���a����	�P�ʁF(m)
	//			RA[ ]		���C���[���a
	//========================================================

	for(i = 0; i < NSEG; ++i){
		RA[i] = wire_radius;	
	}

	//========================================================
	//			�e�G�������g(���C���[)���̃Z�O�����g��
	//			SEGN[ ]
	//========================================================

	////���˔\////
	//���[�v��
	for(hensuu = 1; hensuu <= soshi_no; ++hensuu){
		SEGN[2*hensuu-2] = nbun_lps_o + nbun_l + nbun_lps;
		SEGN[2*hensuu-1] = nbun_lms_o + nbun_l + nbun_lms;
	}
	//���[�v���̐�������
	for(hensuu = 1; hensuu <= soshi_no; ++hensuu){
		SEGN[2*soshi_no + hensuu -1] = nbun_l_o;
	}
	//��������Lx
	for(hensuu = 1; hensuu <= soshi_no-1; ++hensuu){
		SEGN[2*soshi_no + soshi_no + hensuu -1] = nbun_Lx_o;
	}

	////���˔�////
	//���[�v��
	for(hensuu = 1; hensuu <= soshi_no; ++hensuu){
		SEGN[(2*soshi_no + soshi_no + soshi_no -1) + 2*hensuu-2] = nbun_lps_o + nbun_l + nbun_lps;
		SEGN[(2*soshi_no + soshi_no + soshi_no -1) + 2*hensuu-1] = nbun_lms_o + nbun_l + nbun_lms;
	}
	//���[�v���̐�������
	for(hensuu = 1; hensuu <= soshi_no; ++hensuu){
		SEGN[(2*soshi_no + soshi_no + soshi_no -1) + 2*soshi_no + hensuu -1] = nbun_l_o;
	}
	//��������Lx
	for(hensuu = 1; hensuu <= soshi_no-1; ++hensuu){
		SEGN[(2*soshi_no + soshi_no + soshi_no -1) + 2*soshi_no + soshi_no + hensuu -1] = nbun_Lx_o;
	}

	///���d��///
		SEGN[(2*soshi_no + soshi_no + soshi_no -1) + 2*soshi_no + soshi_no + soshi_no -1]
		= (nbun_F_1_o + nbun_h)*2;

	//========================================================
	//			���d�֌W�ݒ�
	//			���d�ʒu		FEDP[ ]
	//			���d�d��		FEDV[ ]
	//
	//			���ʒu  FEDP[ ]
	//					[���d�������ʒu]-[���d���������C���[(���Ԗ�)]
	//			���d��  FEDV[ ]
	//				�@��P���d�̓d��
	//					�s���t���d��p����ꍇ   -2�~(1+j0)[V]
	//					���t���d��p����ꍇ     -1�~(1+j0)[V]
	//				�A���˔t���_�C�|�[����(���˔t�����t���d�A���e�i)
	//					��P���d(���˔\��)   -1�~(1+j0)[V]
	//					��Q���d(���˔���)    1�~(1+j0)[V]
	//						(���˔̕\�Ɨ��ŋ��d�ʑ����t�ɂ���)
	//				�B��P���d���C���[�W�@�̏ꍇ
	//					���̋��d�_���C���[�W�@��K�p����D
	//				�C��P���d���C���[�W�@��p���Ȃ��ꍇ
	//					���̋��d�_���C���[�W�@��p���Ȃ��D
	//				�D�������d�̓d���E�ʑ�������
	//					��P���d�͏����@�A��K�p
	//					��P���d�ȊO�Ő�����s���D
	//========================================================
	FEDP[0] = NSEG - NWIR - (nbun_F_1_o*1 + nbun_h*1) /*+31*//*+31�͉E�̋��d�����狋�d�������Ƃ�*/;				// [�A���e�i�̒���]-[1�{��]
	FEDV[0] = -2.0 * Complex( 1.0 , 0.0 );	// 1.0+j0.0[V]

	/*for(hensuu = 0; hensuu < soshi_no; ++hensuu){
		FEDP[hensuu] = 1;
		FEDV[hensuu] = -1.0 * Complex( 1.0-hensuu*0.1 , 0.0 );
	}
	for(hensuu = 0; hensuu < soshi_no; ++hensuu){
		FEDP[soshi_no+hensuu] = 42*hensuu+426 - hensuu*2-24 -1;
		FEDV[soshi_no+hensuu] = 1.0 * Complex( 1.0-hensuu*0.1 , 0.0 );
	}*/
//��8�f�q�ȊO�ɂ͑Ή����ĂȂ��Ǝv��
	//========================================================
	//			�C���s�[�_���X���׊֌W�ݒ�
	//			���׈ʒu		LOADP[ ]
	//			���גl�@		LOADZ[ ]
	//
	//			���ʒu    LOADP[ ]
	//					[���ׂ������ʒu]-[���ׂ��������C���[(���Ԗ�)]
	//			�����גl�@LOADZ[ ]  (50����ڑ�����ꍇ)
	//				�@���˔ɒ�R��t����ꍇ		 -2�~(50+j0)[��]
	//				�A���C���̓r���ɒ�R��t����ꍇ -1�~(50+j0)[��]
	//				�B�ڑ��Ώۂɂ�炸�C�l�́u�}�C�i�X(-)�v�ɂȂ�
	//								(��R�l�̓x�N�g���ł͖�������)
	//========================================================


	//========================================================
	//			�e�Z�O�����g����̕��ˏ�Ԃ�ݒ�(�������ːݒ�)
	//			RAD_MODE[ ] = 1 or 0
	//				(0 : ���˃J�b�g   1 : �ʏ�̕���)
	//			�������őS��RAD_MODE[ ]=1 �ɐݒ�ς�
	//========================================================
//�S�Ă���˃J�b�g
	for(i = 0; i < ( (nbun_polec+1+nbun_polec1+1+nbun_polec+1+nbun_polec1+1)*soshi_no
			           +(1+nbun_poled+1)*soshi_no
			           +(1+nbun_poled+1)*soshi_no
			           +(1+nbun_Lx+1)*(soshi_no-1) )*2
					   + 1 + nbun_F_1 + nbun_h*2+nbun_F_1+1; ++i){
		RAD_MODE[i] = 0;
	}
	//���˔\
	for(i = (nbun_polec+1+nbun_polec1+1+nbun_polec+1+nbun_polec1+1)*soshi_no
				+(1+nbun_poled+1)*(hensuu); i < (nbun_polec+1+nbun_polec1+1+nbun_polec+1+nbun_polec1+1)*soshi_no
			+(1+nbun_poled+1)*soshi_no
			+(1+nbun_poled+1)*(hensuu); ++i){
		RAD_MODE[i] = 1;
	}

	//���˔�
		for(i = (nbun_polec+1+nbun_polec1+1+nbun_polec+1+nbun_polec1+1)*soshi_no
				+(1+nbun_poled+1)*(hensuu); i < (nbun_polec+1+nbun_polec1+1+nbun_polec+1+nbun_polec1+1)*soshi_no
			+(1+nbun_poled+1)*soshi_no
			+(1+nbun_poled+1)*(hensuu); ++i){
		RAD_MODE[i+(nbun_polec+1+nbun_polec1+1+nbun_polec+1+nbun_polec1+1)*soshi_no
			           +(1+nbun_poled+1)*soshi_no
			           +(1+nbun_poled+1)*soshi_no
			           +(1+nbun_Lx+1)*(soshi_no-1)] = 1;
	}
	//========================================================
	//			�d���E���ˊE�v�Z
	//========================================================	
	ConfCheck();				//�`�󌟍�
	MakeZmn();					//Zmn�쐬
	MakeCurrent();				//�d���v�Z
	MakePhase();				//�d���ʑ��v�Z
	Radiation();				//���ˊE�v�Z
	LapsedTime();				//�o�ߎ��Ԏ擾

	//========================================================
	//			�f�[�^�o��
	//========================================================	
	OutputConf();				//�`��o��
	OutputCurrent();			//�d���o��
	OutputRad();				//���ˊE�o��
	OutputChara();				//�����f�[�^�o��
								//���s

	//========================================================
	//			�C�Ӄf�[�^�o��	free.dat
	//				�g�p�t�@�C���|�C���^ [fp_free]
	//========================================================	
	//fprintf(fp_free," PARA1 = %9.5f    PARA2 = %9.5f\n",PARA1,PARA2);	//�p�����^�[
	//for(i = 0; i < NFED ;++i){
	//	fprintf(fp_free," Zin[%d] = %9.2f %9.2f\n",i,Zin[i].re,Zin[i].im);		//�C���s�[�_���X
	//}
	//fprintf(fp_free," �v�Z���� = %d:%d\n",la_min,la_sec);				//�o�ߎ��ԕ\��
	//fprintf(fp_free," \n\n\n");											//���s
/*
	int nn = 0;
	int zang;
	for(double ii = DegStart ;ii <= DegStart+DegWidth ;ii = ii + DegDelta){
		if( (int)(ii*10000000)==-4.0*10000000){
			zang = nn;
		}
		nn = nn + 1;
	}

	fprintf(fp_free," %9.3f %7.3f **\n",RLGAI[zang],AR[zang]);*/
	double F_2_pole=(0.125+0.001*(-25)) * RAMDA0;
	fprintf(fp_free,"%d",int(F_2_pole/segtyou+0.5));
	//========================================================
	//			�z��̊J��
	//========================================================	
	DelMatAll();				//�z��J��(���R�����g�A�E�g�֎~��)
	}
}
}
}
}}
//============================================================================
//					�ϐ��Ɣz��ꗗ
//============================================================================
//
//--�`�� ���d
//		double	FREQ0					�݌v���g��f0
//		double	RAMDA0					�݌v���g���̎��R��Ԕg��
//		double	USEF					���d���g��
//		int		NWIR					�S���C���[��
//		int		NSEG					�S�Z�O�����g��
//		int		NFED					�S���d�_��
//		int		NSEG0					�S�d���v�Z�_��
//		int		NLOAD					��R�����א�
//--���ˊE�v�Z
//		int		AXMODE   				(0=�ӌŒ�, 1=�ƌŒ�)	
//		double	DegDelta 				���ݕ�
//		double	DegStart 				�����p
//		double	DegWidth 				�͈�
//		double	FixAngle 				�Œ莲�p�x
//--���ԑ���
//		int		la_min	la_sec			�� �b
//--�v�Z�񂵗p�p�����[�^
//		double	PARA1,PARA2				�p�����[�^�@�ƇA
//--�`��z��
//		double	RX[ ]  RY[ ]  RZ[ ]		�e�Z�O�����g�̎n�_
//		double�@SX[ ]  SY[ ]  SZ[ ]		�P�ʃx�N�g���̐���
//		double	SEGL[ ]	RA[ ]			�e�Z�O�����g�̒��� ���a�@		
//		int�@	SEGN[ ]					�e���C���[�̃Z�O�����g��
//		int		FEDP[ ]					�e���d�ʒu
//		complex	FEDV[ ]					�e���d�d��
//		int		RAD_MODE[ ]				�e�Z�O�����g�̕��˃��[�h
//		int		LOADP[ ]				�e�C���s�[�_���X���׈ʒu
//		complex LOADZ[ ]				�e�C���s�[�_���X���גl
//--�d���z��
//		complex Zmn[ ][ ]				�C���s�[�_���X�s��	
//		complex Im[ ]					�d�����z
//		double  PhaseIm[ ]				�d���̈ʑ�
//--���̓C���s�[�_���X
//		complex Zin[ ]					�e���d�_�̓��̓C���s�[�_���X
//		double  VSWR_50[ ]  VSWR_75[ ]	�e���d�_��50����75���ɑ΂���VSWR
//--���ˊE�z��
//		complex TRAD[ ] FRAD[ ]			E�Ƃ�E�ӂ̕��ˊE
//		complex RRAD[ ] LRAD[ ]			ER ��EL �̕��ˊE
//		double TPHASE[ ] FPHASE[ ]		E�Ƃ�E�ӂ̈ʑ�
//		double RPHASE[ ] LPHASE[ ]		ER ��EL �̈ʑ�
//		double TGAI[ ] FGAI[ ]			E�Ƃ�E�ӂ̗���
//		double TFGAI[ ]					E�Ƃ�E�ӂ̍��v�̗���
//		double RLGAI[ ]					ER ��EL �̍��v�̗���
//		double RGAI[ ] LGAI[ ]			ER ��EL �̗���
//		double TPAT[ ] FPAT[ ]			E�Ƃ�E�ӂ̃p�^�[��
//		double RPAT[ ] LPAT[ ]			ER ��EL �̃p�^�[��
//		double AR[ ]					����
//--��{�����萔
//		double PI						��
//		double C						����
//		double e0						�^�󒆂̓�������0
//		double u0						�^�󒆂̗U�d����0
//		double R						1.0 �Œ�
//--���ȗ��p
//		complex J						�����P�ʁ@J=��-1
//--�ϕ��p
//		int GaussTenNo					�K�E�X�ϕ����_  [ 4�Œ�]�i�ʏ�p�j
//		int GaussTenSpe					�K�E�X�ϕ����_��[40�Œ�]�i���ٓ_�p�j
//--�v�Z�p
//		double k0						k0
//		complex uair					u(air)
//		double Beta;					��
//============================================================================
//============================================================================
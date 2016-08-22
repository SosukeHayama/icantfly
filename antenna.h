//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
//			アンテナ形状ヘッダー　　antenna.h
//			貫通形アンテナ
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
//				計算
//////////////////////////////////////////////////////////////////////////////////
void Calculation(void)
{
	for(int AA = 59;AA<= 59; ++AA){		//x軸からの距離d用(固定するときver)
										//d_pole = 0.001 * RAMDA0 * AA;
	for(int CC = 0; CC<= 0; ++CC){		//f特用
	for(int DD = 125; DD<= 125; ++DD){	//F～#1間の直線導体長最適化用
										//F_1_pole = 0.001*DD * RAMDA0;1.169
	for(int EE = 1071; EE<= 1071; ++EE){	//C最適化用
	for(int BB = 1023; BB<= 1023; ++BB){	//Lx+C/4最適化用
										//Lx_pole = 0.001*BB * RAMDA0;
	//========================================================
	//			ローカル変数
	//========================================================
	//--カウンタ
	int i,j;

	//--計算で使用する変数
	double rrx , rry , rrz;				//セグメント始点
	double drx , dry , drz;				//セグメント終点
	double ssx , ssy , ssz;				//セグメントの単位ベクトル
	double dcl;							//セグメント長さ

	//--エレメント数関係（ワイヤー本数）
	int n_array;

	//--各部分の分割数
	int nbun_l , nbun_l_o;			//正方形の1辺の分割数
	int nbun_d;						//x軸と水平導体の間隔d分の分割数
	int nbun_lps , nbun_lps_o;		//(正方形の1辺の長さ/2 + 間隔d)の分割数
	int nbun_lms , nbun_lms_o;		//(正方形の1辺の長さ/2 - 間隔d)の分割数
	int nbun_h , nbun_h_o;			//給電線垂直部の分割数
	int nbun_F_1 , nbun_F_1_o;		//ループ#1と給電線の間のライン
	int nbun_Lx , nbun_Lx_o;		//直線導体Lx
	int nbun_F;						//給電線水平部縦


	//--各部のセグメント数
	int n_pole;

	//--設計値関係
	double wire_radius;			//ワイヤー半径
	double C_pole;				//正方形ループ周囲長
	double l_pole;				//正方形の1辺の長さ
	double d_pole;				//x軸と水平導体の間隔d
	double lpd_pole;			//正方形の1辺の長さ/2 + 間隔d
	double lmd_pole;			//正方形の1辺の長さ/2 - 間隔d
	double h_pole;				//高さh
	double F_1_pole;			//ループ#1と給電線の間のライン
	double Lx_pole;				//直線導体Lx
	double F_pole;				//給電線水平部縦

	int hensuu;					//なんにでも使っていい変数
	int hensuu2;

	int soshi_no = 4;

	//========================================================
	//			放射界計算設定		単位：(deg)
	//========================================================
	//--計算モード
	AXMODE = 0;			//(0=φ固定, 1=θ固定)

	//--可変軸
	DegDelta =    1.0;	//刻み幅("0.0"禁止)
	DegStart =  -90.0;	//初期角
	DegWidth =  180.0;	//範囲

	//--固定軸
	FixAngle =  0.0;	//固定角度

	//========================================================
	//			設計周波数		FREQ0	(f0)	単位：(Hz)
	//			給電周波数		USEF	(f)		単位：(Hz)
	//========================================================
	//--設計周波数設定　
	FREQ0 = 3.081 * pow(10.0,9.0);

	//--電源周数設定
	//USEF  = FREQ0 * 1.0;					//f特以外のとき
	USEF  = FREQ0 * (1.000+0.001*CC);					//f特のとき

	//--設計周波数の自由空間波長
	RAMDA0= C/FREQ0;

	//========================================================
	//			形状値定義
	//			長さはRAMDA0を基準　係数で値定義(○○λ)
	//========================================================

	C_pole = 0.001 * EE * RAMDA0;				//正方形ループ周囲長
	l_pole = C_pole/4.0;				//正方形の1辺の長さ

	d_pole = 0.001 * RAMDA0 * AA;			//x軸と水平導体の間隔d
	lpd_pole = l_pole/2.0 + d_pole;		//正方形の1辺の長さ/2 + 間隔d
	lmd_pole = l_pole/2.0 - d_pole;		//正方形の1辺の長さ/2 - 間隔d

	h_pole = 0.125 * RAMDA0;			//高さh

	F_1_pole = 0.001*DD * RAMDA0;			//ループ#1と給電線の間のライン
	Lx_pole = 0.001*BB * RAMDA0 - C_pole/4.0;			//直線導体Lx
	F_pole = 0.125 * RAMDA0;			//給電線水平部縦

	wire_radius = 0.005 * RAMDA0;	//ワイヤー半径

	//========================================================
	//			出力ファイル用のパラメータ
	//			PARA1 と PARA2 に値を入力する．
	//========================================================
	PARA1 = d_pole/RAMDA0;
	PARA2 = C_pole / RAMDA0;
	PARA3 = (Lx_pole+l_pole)/RAMDA0;
	
	//PARA1 = CC;
	PARA4 = USEF/FREQ0;


	//========================================================
	//			全ワイヤー数	NWIR
	//			セグメント数	NSEG
	//			給電点数		NFED
	//			全電流計算点数	NSEG0
	//			Ｚ装荷数		NLOAD		
	//========================================================

	//--全ワイヤー数設定
	n_array = (2*soshi_no					//ループ
			+ 1*soshi_no					//水平導体
			+ 1*(soshi_no-1))*2				//直線導体Lx
			+ 1;							//給電線
	NWIR = n_array;
  	double segtyou=0.025*RAMDA0;
	
	//--分割数入力
	nbun_l = int(l_pole/segtyou+0.5);							//正方形の1辺の分割数
	nbun_d = int(d_pole/segtyou+0.5);								//x軸と水平導体の間隔d分の分割数
	nbun_lps = int(lpd_pole/segtyou+0.5);			//(正方形の1辺の長さ/2 + 間隔d)の分割数
	nbun_lms = int(lmd_pole/segtyou+0.5);			//(正方形の1辺の長さ/2 - 間隔d)の分割数
	nbun_h = int(h_pole/segtyou+0.5);								//垂直導線F-F'の分割数
	nbun_F_1 = int(F_1_pole/segtyou+0.5);							//ループ#1と給電線の間のライン
	nbun_Lx = int(Lx_pole/segtyou+0.5);							//直線導体Lx
	nbun_F = int(F_pole/segtyou+0.5);								//給電線水平部縦

	//--オーバーラップ込分割数
	nbun_l_o = nbun_l + 2;					//正方形の1辺の分割数
	nbun_lps_o = nbun_lps + 1;				//(正方形の1辺の長さ/2 + 間隔d)の分割数
	nbun_lms_o = nbun_lms + 1;				//(正方形の1辺の長さ/2 - 間隔d)の分割数
	nbun_h_o = nbun_h + 1;					//垂直導線F-F'の分割数
	nbun_F_1_o = nbun_F_1 + 1;				//ループ#1と給電線の間のライン
	nbun_Lx_o = nbun_Lx + 2;				//直線導体Lx


	//--全セグメント数設定
	n_pole = ( (nbun_lps_o + nbun_l + nbun_lps + nbun_lms_o + nbun_l + nbun_lms)*soshi_no	//ループ
			+ nbun_l_o*soshi_no																//水平導体
			+ nbun_Lx_o*(soshi_no-1) )*2													//直線導体Lx
			+ (nbun_F_1_o + nbun_h)*2;

	NSEG = n_pole;						//全セグメント数

	//--全電流計算点数  [NSEG0=NSEG-NWIR]	※書き換え＆削除　不可
	NSEG0 = NSEG - NWIR;

	//--給電点数
	NFED = 1;

	//--インピーダンス装荷数
	NLOAD = 0;

	//========================================================
	//			配列確保と初期化
	//			[放射界計算設定,NWIR,NSEG,NFED,NLOAD]
	//				設定後に関数呼び出し
	//			※書き換え＆削除　不可
	//========================================================
	MakeMatAll();								//配列確保(※コメントアウト禁止※)
	Initialization();							//配列と計算用変数の初期化
	
	//========================================================
	//			設計値入力		単位：(m)
	//			セグメント始点	RX[ ] , RY[ ] , RZ[ ] 
	//			単位ベクトル	SX[ ] , SY[ ] , SZ[ ]
	//			セグメント長	SEGL[ ]
	//========================================================

	////反射板表////
	//反射板_表_ループ_上部_オーバーラップ
		i = 0;
		//--始点 終点 セグメント長計算
		rrx = F_1_pole;							//始点位置計算x
		rry = -(lmd_pole/nbun_lms) - d_pole;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = F_1_pole;							//始点位置計算x
		dry = -d_pole;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長

	//反射板_表_ループ_上部_左辺
	for(j = 0; j < nbun_lps_o - 1; ++j){
		i = j + (1);
		//--始点 終点 セグメント長計算
		rrx = F_1_pole;							//始点位置計算x
		rry = (lpd_pole/nbun_lps)*(j+0) - d_pole;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = F_1_pole;							//始点位置計算x
		dry = (lpd_pole/nbun_lps)*(j+1) - d_pole;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長
	}

	//反射板_表_ループ_上部_上辺
	for(j = 0; j < nbun_l; ++j){
		i = j + (nbun_lps_o);
		//--始点 終点 セグメント長計算
		rrx = (l_pole/nbun_l)*(j+0) + F_1_pole;							//始点位置計算x
		rry = l_pole/2.0;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = (l_pole/nbun_l)*(j+1) + F_1_pole;							//始点位置計算x
		dry = l_pole/2.0;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長
	}

	//反射板_表_ループ_上部_右辺
	for(j = 0; j < nbun_lps; ++j){
		i = j + (nbun_lps_o + nbun_l);
		//--始点 終点 セグメント長計算
		rrx = F_1_pole + l_pole;							//始点位置計算x
		rry = -(lpd_pole/nbun_lps)*(j+0) + l_pole/2.0;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = F_1_pole + l_pole;							//始点位置計算x
		dry = -(lpd_pole/nbun_lps)*(j+1) + l_pole/2.0;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長
	}

	//反射板_表_ループ_下部_オーバーラップ
		i = (nbun_lps_o + nbun_l + nbun_lps);
		//--始点 終点 セグメント長計算
		rrx = F_1_pole + l_pole;							//始点位置計算x
		rry = (lpd_pole/nbun_lps) - d_pole;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = F_1_pole + l_pole;							//始点位置計算x
		dry = -d_pole;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長

	//反射板_表_ループ_下部_右辺
	for(j = 0; j < nbun_lms_o - 1; ++j){
		i = j + (nbun_lps_o + nbun_l + nbun_lps) + (1);
		//--始点 終点 セグメント長計算
		rrx = F_1_pole + l_pole;							//始点位置計算x
		rry = - (lmd_pole/nbun_lms)*(j+0) -d_pole;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = F_1_pole + l_pole;							//始点位置計算x
		dry = - (lmd_pole/nbun_lms)*(j+1) -d_pole;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長
	}
	
	//反射板_表_ループ_下部_下辺
	for(j = 0; j < nbun_l; ++j){
		i = j + (nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o);
		//--始点 終点 セグメント長計算
		rrx = -(l_pole/nbun_l)*(j+0) + F_1_pole + l_pole;							//始点位置計算x
		rry = -l_pole/2.0;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = -(l_pole/nbun_l)*(j+1) + F_1_pole + l_pole;							//始点位置計算x
		dry = -l_pole/2.0;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長
	}

	//反射板_表_ループ_下部_左辺
	for(j = 0; j < nbun_lms; ++j){
		i = j + (nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l);
		//--始点 終点 セグメント長計算
		rrx = F_1_pole;							//始点位置計算x
		rry = (lmd_pole/nbun_lms)*(j+0) -l_pole/2.0;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = F_1_pole;							//始点位置計算x
		dry = (lmd_pole/nbun_lms)*(j+1) -l_pole/2.0;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長
	}

	//反射板_表_ループ#N
	for(hensuu = 2; hensuu <= soshi_no; ++hensuu){
		for(j = 0; j < (nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms); ++j){
			i = j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*(hensuu-1);
			//--計算結果代入
			RX[i] = RX[j] + (l_pole + Lx_pole)*(hensuu-1);				//始点位置x
			RY[i] = RY[j];											//始点位置y
			RZ[i] = RZ[j];											//始点位置z
			SX[i] = SX[j];											//単位ベクトルx
			SY[i] = SY[j];											//単位ベクトルy
			SZ[i] = SZ[j];											//単位ベクトルz
			SEGL[i] = SEGL[j];										//セグメント長
		}
	}

	//反射板_表_水平導線_オーバーラップ
		i = ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no;
		//--始点 終点 セグメント長計算
		rrx = F_1_pole;							//始点位置計算x
		rry = (lpd_pole/nbun_lps)-d_pole;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = F_1_pole;							//始点位置計算x
		dry = -d_pole;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長

	//反射板_表_水平導線
	for(j = 0; j < nbun_l; ++j){
		i = j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ 1;
		//--始点 終点 セグメント長計算
		rrx = (l_pole/nbun_l)*(j+0) + F_1_pole;							//始点位置計算x
		rry = -d_pole;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx =(l_pole/nbun_l)*(j+1) + F_1_pole;								//始点位置計算x
		dry = -d_pole;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長
	}

	//反射板_表_水平導線_オーバーラップ
		i = ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l + 1;
		//--始点 終点 セグメント長計算
		rrx = F_1_pole + l_pole;							//始点位置計算x
		rry = -d_pole;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = F_1_pole + l_pole;							//始点位置計算x
		dry = -(lmd_pole/nbun_lms)-d_pole;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長

	//反射板_表_水平導体_#2～#N
	for(hensuu = 2; hensuu <= soshi_no; ++hensuu){
		for(j = 0; j < nbun_l_o; ++j){
			i = j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
				+ nbun_l_o*(hensuu-1);
			//--計算結果代入
			RX[i] = RX[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no] + (l_pole + Lx_pole)*(hensuu-1);				//始点位置x
			RY[i] = RY[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no];											//始点位置y
			RZ[i] = RZ[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no];											//始点位置z
			SX[i] = SX[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no];											//単位ベクトルx
			SY[i] = SY[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no];											//単位ベクトルy
			SZ[i] = SZ[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no];											//単位ベクトルz
			SEGL[i] = SEGL[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no];										//セグメント長
		}
	}

	//反射板_表_直線導体Lx_オーバーラップ
		i = ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no;
		//--始点 終点 セグメント長計算
		rrx = -(l_pole/nbun_l) + F_1_pole + l_pole;							//始点位置計算x
		rry = -d_pole;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = F_1_pole + l_pole;							//始点位置計算x
		dry = -d_pole;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長

	//反射板_表_直線導体Lx
	for(j = 0; j < nbun_Lx; ++j){
		i = j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ 1;
		//--始点 終点 セグメント長計算
		rrx = (Lx_pole/nbun_Lx)*(j+0) + F_1_pole + l_pole;							//始点位置計算x
		rry = -d_pole;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx =(Lx_pole/nbun_Lx)*(j+1) + F_1_pole + l_pole;								//始点位置計算x
		dry = -d_pole;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長
	}

	//反射板_表_直線導体Lx_オーバーラップ
		i = ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx + 1;
		//--始点 終点 セグメント長計算
		rrx = F_1_pole + l_pole + Lx_pole;							//始点位置計算x
		rry = -d_pole;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = (l_pole/nbun_l) + F_1_pole + l_pole + Lx_pole;							//始点位置計算x
		dry = -d_pole;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長

	//反射板_表_直線導体Lx_残り
	for(hensuu = 2; hensuu <= soshi_no-1; ++hensuu){
		for(j = 0; j < nbun_Lx_o; ++j){
		i = j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx_o*(hensuu-1);
			//--計算結果代入
			RX[i] = RX[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no] + (l_pole + Lx_pole)*(hensuu-1);				//始点位置x
			RY[i] = RY[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no];											//始点位置y
			RZ[i] = RZ[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no];											//始点位置z
			SX[i] = SX[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no];											//単位ベクトルx
			SY[i] = SY[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no];											//単位ベクトルy
			SZ[i] = SZ[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no];											//単位ベクトルz
			SEGL[i] = SEGL[j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no + nbun_l_o*soshi_no];										//セグメント長
		}
	}

	//反射板裏
		for(j = 0; j < ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
						+ nbun_l_o*soshi_no
						+ nbun_Lx_o*(soshi_no-1); ++j){
			i = j + ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
				+ nbun_l_o*soshi_no
				+ nbun_Lx_o*(soshi_no-1);
			//--計算結果代入
			RX[i] = RX[j];											//始点位置x
			RY[i] = RY[j];											//始点位置y
			RZ[i] = -RZ[j];											//始点位置z
			SX[i] = SX[j];											//単位ベクトルx
			SY[i] = SY[j];											//単位ベクトルy
			SZ[i] = SZ[j];											//単位ベクトルz
			SEGL[i] = SEGL[j];										//セグメント長
		}

	//給電線_オーバーラップ
		i = ( ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx_o*(soshi_no-1) )*2;
		//--始点 終点 セグメント長計算
		rrx = (l_pole/nbun_l) + F_1_pole;							//始点位置計算x
		rry = -d_pole;							//始点位置計算y
		rrz = -h_pole;														//始点位置計算z
		drx = F_1_pole;							//始点位置計算x
		dry = -d_pole;							//始点位置計算y
		drz = -h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長

	//給電線_F～#1_反射板裏
	for(j = 0; j < nbun_F_1; ++j){
		i = j + ( ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx_o*(soshi_no-1) )*2
			+ 1;
		//--始点 終点 セグメント長計算
		rrx = -(F_1_pole/nbun_F_1)*(j+0) + F_1_pole;							//始点位置計算x
		rry = -d_pole;							//始点位置計算y
		rrz = -h_pole;														//始点位置計算z
		drx = -(F_1_pole/nbun_F_1)*(j+1) + F_1_pole;								//始点位置計算x
		dry = -d_pole;							//始点位置計算y
		drz = -h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長
	}

	//給電線_垂直部
	for(j = 0; j < nbun_h*2; ++j){
		i = j + ( ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx_o*(soshi_no-1) )*2
			+ nbun_F_1_o;
		//--始点 終点 セグメント長計算
		rrx = 0;							//始点位置計算x
		rry = - d_pole;							//始点位置計算y
		rrz = (h_pole/nbun_h)*(j+0)-h_pole;														//始点位置計算z
		drx = 0;								//始点位置計算x
		dry = - d_pole;							//始点位置計算y
		drz = (h_pole/nbun_h)*(j+1)-h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長
	}

	//給電線_F～#1_反射板表
	for(j = 0; j < nbun_F_1; ++j){
		i = j + ( ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx_o*(soshi_no-1) )*2
			+ nbun_F_1_o + nbun_h*2;
		//--始点 終点 セグメント長計算
		rrx = (F_1_pole/nbun_F_1)*(j+0);							//始点位置計算x
		rry = -d_pole;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = (F_1_pole/nbun_F_1)*(j+1);								//始点位置計算x
		dry = -d_pole;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長
	}

	//給電線_オーバーラップ
		i = ( ((nbun_lps_o + nbun_l + nbun_lps) + (nbun_lms_o + nbun_l + nbun_lms))*soshi_no
			+ nbun_l_o*soshi_no
			+ nbun_Lx_o*(soshi_no-1) )*2
			+ nbun_F_1_o + nbun_h*2 + nbun_F_1;
		//--始点 終点 セグメント長計算
		rrx = F_1_pole;							//始点位置計算x
		rry = -d_pole;							//始点位置計算y
		rrz = h_pole;														//始点位置計算z
		drx = (l_pole/nbun_l) + F_1_pole;							//始点位置計算x
		dry = -d_pole;							//始点位置計算y
		drz = h_pole;														//終点位置計算z
		//--セグメント長 単位方向ベクトルの各成分の計算
		dcl = sqrt( pow( (rrx-drx),2 )
				  + pow( (rry-dry),2 )
				  + pow( (rrz-drz),2 ) );						//セグメント長計算
		ssx = (drx-rrx) / dcl;									//単位ベクトルx
		ssy = (dry-rry) / dcl;									//単位ベクトルy
		ssz = (drz-rrz) / dcl;									//単位ベクトルz
		//--計算結果代入
		RX[i] = rrx;											//始点位置x
		RY[i] = rry;											//始点位置y
		RZ[i] = rrz;											//始点位置z
		SX[i] = ssx;											//単位ベクトルx
		SY[i] = ssy;											//単位ベクトルy
		SZ[i] = ssz;											//単位ベクトルz
		SEGL[i] = dcl;											//セグメント長
	
	//========================================================
	//			各セグメントのワイヤー半径入力	単位：(m)
	//			RA[ ]		ワイヤー半径
	//========================================================

	for(i = 0; i < NSEG; ++i){
		RA[i] = wire_radius;	
	}

	//========================================================
	//			各エレメント(ワイヤー)毎のセグメント数
	//			SEGN[ ]
	//========================================================

	////反射板表////
	//ループ部
	for(hensuu = 1; hensuu <= soshi_no; ++hensuu){
		SEGN[2*hensuu-2] = nbun_lps_o + nbun_l + nbun_lps;
		SEGN[2*hensuu-1] = nbun_lms_o + nbun_l + nbun_lms;
	}
	//ループ内の水平導体
	for(hensuu = 1; hensuu <= soshi_no; ++hensuu){
		SEGN[2*soshi_no + hensuu -1] = nbun_l_o;
	}
	//垂直導体Lx
	for(hensuu = 1; hensuu <= soshi_no-1; ++hensuu){
		SEGN[2*soshi_no + soshi_no + hensuu -1] = nbun_Lx_o;
	}

	////反射板裏////
	//ループ部
	for(hensuu = 1; hensuu <= soshi_no; ++hensuu){
		SEGN[(2*soshi_no + soshi_no + soshi_no -1) + 2*hensuu-2] = nbun_lps_o + nbun_l + nbun_lps;
		SEGN[(2*soshi_no + soshi_no + soshi_no -1) + 2*hensuu-1] = nbun_lms_o + nbun_l + nbun_lms;
	}
	//ループ内の水平導体
	for(hensuu = 1; hensuu <= soshi_no; ++hensuu){
		SEGN[(2*soshi_no + soshi_no + soshi_no -1) + 2*soshi_no + hensuu -1] = nbun_l_o;
	}
	//垂直導体Lx
	for(hensuu = 1; hensuu <= soshi_no-1; ++hensuu){
		SEGN[(2*soshi_no + soshi_no + soshi_no -1) + 2*soshi_no + soshi_no + hensuu -1] = nbun_Lx_o;
	}

	///給電線///
		SEGN[(2*soshi_no + soshi_no + soshi_no -1) + 2*soshi_no + soshi_no + soshi_no -1]
		= (nbun_F_1_o + nbun_h)*2;

	//========================================================
	//			給電関係設定
	//			給電位置		FEDP[ ]
	//			給電電圧		FEDV[ ]
	//
	//			※位置  FEDP[ ]
	//					[給電したい位置]-[給電したいワイヤー(何番目)]
	//			※電圧  FEDV[ ]
	//				①第１給電の電圧
	//					不平衡給電を用いる場合   -2×(1+j0)[V]
	//					平衡給電を用いる場合     -1×(1+j0)[V]
	//				②反射板付きダイポール等(反射板付き平衡給電アンテナ)
	//					第１給電(反射板表側)   -1×(1+j0)[V]
	//					第２給電(反射板裏側)    1×(1+j0)[V]
	//						(反射板の表と裏で給電位相を逆にする)
	//				③第１給電がイメージ法の場合
	//					他の給電点もイメージ法を適用する．
	//				④第１給電がイメージ法を用いない場合
	//					他の給電点もイメージ法を用いない．
	//				⑤複数給電の電圧・位相差制御
	//					第１給電は条件①②を適用
	//					第１給電以外で制御を行う．
	//========================================================
	FEDP[0] = NSEG - NWIR - (nbun_F_1_o*1 + nbun_h*1) /*+31*//*+31は右の給電線から給電したいとき*/;				// [アンテナの中央]-[1本目]
	FEDV[0] = -2.0 * Complex( 1.0 , 0.0 );	// 1.0+j0.0[V]

	/*for(hensuu = 0; hensuu < soshi_no; ++hensuu){
		FEDP[hensuu] = 1;
		FEDV[hensuu] = -1.0 * Complex( 1.0-hensuu*0.1 , 0.0 );
	}
	for(hensuu = 0; hensuu < soshi_no; ++hensuu){
		FEDP[soshi_no+hensuu] = 42*hensuu+426 - hensuu*2-24 -1;
		FEDV[soshi_no+hensuu] = 1.0 * Complex( 1.0-hensuu*0.1 , 0.0 );
	}*/
//↑8素子以外には対応してないと思う
	//========================================================
	//			インピーダンス装荷関係設定
	//			装荷位置		LOADP[ ]
	//			装荷値　		LOADZ[ ]
	//
	//			※位置    LOADP[ ]
	//					[装荷したい位置]-[装荷したいワイヤー(何番目)]
	//			※装荷値　LOADZ[ ]  (50Ωを接続する場合)
	//				①反射板に抵抗を付ける場合		 -2×(50+j0)[Ω]
	//				②ラインの途中に抵抗を付ける場合 -1×(50+j0)[Ω]
	//				③接続対象によらず，値は「マイナス(-)」になる
	//								(抵抗値はベクトルでは無いため)
	//========================================================


	//========================================================
	//			各セグメントからの放射状態を設定(部分放射設定)
	//			RAD_MODE[ ] = 1 or 0
	//				(0 : 放射カット   1 : 通常の放射)
	//			初期化で全てRAD_MODE[ ]=1 に設定済み
	//========================================================
//全てを放射カット
	for(i = 0; i < ( (nbun_polec+1+nbun_polec1+1+nbun_polec+1+nbun_polec1+1)*soshi_no
			           +(1+nbun_poled+1)*soshi_no
			           +(1+nbun_poled+1)*soshi_no
			           +(1+nbun_Lx+1)*(soshi_no-1) )*2
					   + 1 + nbun_F_1 + nbun_h*2+nbun_F_1+1; ++i){
		RAD_MODE[i] = 0;
	}
	//反射板表
	for(i = (nbun_polec+1+nbun_polec1+1+nbun_polec+1+nbun_polec1+1)*soshi_no
				+(1+nbun_poled+1)*(hensuu); i < (nbun_polec+1+nbun_polec1+1+nbun_polec+1+nbun_polec1+1)*soshi_no
			+(1+nbun_poled+1)*soshi_no
			+(1+nbun_poled+1)*(hensuu); ++i){
		RAD_MODE[i] = 1;
	}

	//反射板裏
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
	//			電流・放射界計算
	//========================================================	
	ConfCheck();				//形状検査
	MakeZmn();					//Zmn作成
	MakeCurrent();				//電流計算
	MakePhase();				//電流位相計算
	Radiation();				//放射界計算
	LapsedTime();				//経過時間取得

	//========================================================
	//			データ出力
	//========================================================	
	OutputConf();				//形状出力
	OutputCurrent();			//電流出力
	OutputRad();				//放射界出力
	OutputChara();				//特性データ出力
								//改行

	//========================================================
	//			任意データ出力	free.dat
	//				使用ファイルポインタ [fp_free]
	//========================================================	
	//fprintf(fp_free," PARA1 = %9.5f    PARA2 = %9.5f\n",PARA1,PARA2);	//パラメター
	//for(i = 0; i < NFED ;++i){
	//	fprintf(fp_free," Zin[%d] = %9.2f %9.2f\n",i,Zin[i].re,Zin[i].im);		//インピーダンス
	//}
	//fprintf(fp_free," 計算時間 = %d:%d\n",la_min,la_sec);				//経過時間表示
	//fprintf(fp_free," \n\n\n");											//改行
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
	//			配列の開放
	//========================================================	
	DelMatAll();				//配列開放(※コメントアウト禁止※)
	}
}
}
}
}}
//============================================================================
//					変数と配列一覧
//============================================================================
//
//--形状 給電
//		double	FREQ0					設計周波数f0
//		double	RAMDA0					設計周波数の自由空間波長
//		double	USEF					給電周波数
//		int		NWIR					全ワイヤー数
//		int		NSEG					全セグメント数
//		int		NFED					全給電点数
//		int		NSEG0					全電流計算点数
//		int		NLOAD					抵抗等装荷数
//--放射界計算
//		int		AXMODE   				(0=φ固定, 1=θ固定)	
//		double	DegDelta 				刻み幅
//		double	DegStart 				初期角
//		double	DegWidth 				範囲
//		double	FixAngle 				固定軸角度
//--時間測定
//		int		la_min	la_sec			分 秒
//--計算回し用パラメータ
//		double	PARA1,PARA2				パラメータ①と②
//--形状配列
//		double	RX[ ]  RY[ ]  RZ[ ]		各セグメントの始点
//		double　SX[ ]  SY[ ]  SZ[ ]		単位ベクトルの成分
//		double	SEGL[ ]	RA[ ]			各セグメントの長さ 半径　		
//		int　	SEGN[ ]					各ワイヤーのセグメント数
//		int		FEDP[ ]					各給電位置
//		complex	FEDV[ ]					各給電電圧
//		int		RAD_MODE[ ]				各セグメントの放射モード
//		int		LOADP[ ]				各インピーダンス装荷位置
//		complex LOADZ[ ]				各インピーダンス装荷値
//--電流配列
//		complex Zmn[ ][ ]				インピーダンス行列	
//		complex Im[ ]					電流分布
//		double  PhaseIm[ ]				電流の位相
//--入力インピーダンス
//		complex Zin[ ]					各給電点の入力インピーダンス
//		double  VSWR_50[ ]  VSWR_75[ ]	各給電点の50Ωと75Ωに対するVSWR
//--放射界配列
//		complex TRAD[ ] FRAD[ ]			EθとEφの放射界
//		complex RRAD[ ] LRAD[ ]			ER とEL の放射界
//		double TPHASE[ ] FPHASE[ ]		EθとEφの位相
//		double RPHASE[ ] LPHASE[ ]		ER とEL の位相
//		double TGAI[ ] FGAI[ ]			EθとEφの利得
//		double TFGAI[ ]					EθとEφの合計の利得
//		double RLGAI[ ]					ER とEL の合計の利得
//		double RGAI[ ] LGAI[ ]			ER とEL の利得
//		double TPAT[ ] FPAT[ ]			EθとEφのパターン
//		double RPAT[ ] LPAT[ ]			ER とEL のパターン
//		double AR[ ]					軸比
//--基本物理定数
//		double PI						π
//		double C						光速
//		double e0						真空中の透磁率ε0
//		double u0						真空中の誘電率μ0
//		double R						1.0 固定
//--式簡略用
//		complex J						虚数単位　J=√-1
//--積分用
//		int GaussTenNo					ガウス積分分点  [ 4固定]（通常用）
//		int GaussTenSpe					ガウス積分分点数[40固定]（特異点用）
//--計算用
//		double k0						k0
//		complex uair					u(air)
//		double Beta;					β
//============================================================================
//============================================================================
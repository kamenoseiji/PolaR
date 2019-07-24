########
# To use readSAM45.R, apply source('date.R') in advance.
#
#-------- Math Constants
DEGRAD <- pi/180.0
#-------- Constants
ARYMAX <- 32	# アレイ数の最大値
BINMAX <- 4		# ビン分割の最大値
SPWMAX <- 2		# スペクトラルウインドウ指定の最大値
SUBARY <- 2		# サブアレイ分割の最大値
CHMAX <- 4096	# ロギングチャネル数の最大値

#-------- Read 1st header in SAM45 logging
readSAM45head1 <- function( file_ptr ){
	head1 <- list(
		readChar(file_ptr, nchars=4),			# LC : P+O ロギングヘッダータイプ L0, L1, L2, ED
		readBin(file_ptr, integer(), size=4),	# LC : P+O レコード長
		readChar(file_ptr, nchars=8),			# LC : P+O ヘッダー情報のバージョン V004
		readChar(file_ptr, nchars=64),			# LC : P+O ロギングファイル名 SAM45.指示書名.group.project.ロギングID
		readChar(file_ptr, nchars=8))			# LC : P+O バックエンド種別 SAM45
	names(head1) <- c("crec_type", "irec_len", "cversion", "clog_name", "cbe_type")
	return(head1)
}

#-------- Read 2nd header in SAM45 logging
readSAM45head2 <- function( file_ptr ){
	namelist <- c("crec_type", "irec_len", "dant_off", "dsub_off", "dsrc_pos", "dstrv", "dmap_pos", "dmult_off", "dbsw_freq", "dsw_freq", "dflif", "dcent_freq", "dtrk_freq", "dbebw", "dfqdat_f0", "dfqdat_fq", "dfqdat_ch", "dbechwid", "dberes", "deffa", "deffb", "deffl", "defss", "dgain", "dhpbw", "dpoldr", "dbin_dwell", "diptim", "reserve8", "ich_bind", "ich_range", "icalb_int", "ichannel", "iipint", "ipordr", "ichan_avg", "istart_chan", "iend_chan", "drms_value", "on_axis", "irms_judge", "irms_mode", "irms_onnum", "reserve4", "imd_obs", "icoord_type", "iscan_coord", "isw_mod", "ivdef", "ivref", "imlt_no", "iary_ifatt", "iary_refno", "iary_usefg", "ibin_num", "in_spec_window", "clog_id", "cst_time", "cobs_file", "cobj_name", "csrc_comment", "cobs_user", "cgroup", "cproject", "ctrk_type", "cepoch", "cscan_coord", "cseq_ptn", "csid_type", "cfe_type", "chorn_typ", "cpol_typ", "csub_array", "cinput_mode", "cbut_scale", "crequant_rescale", "crequant_scale", "cmake_hfs_table", "cmisc_flags_hfs", "cmisc_flags_lfs", "cmisc_flags_mult", "cfreq_prof_synth", "cwin_func_id", "crms_array", "i2beam_mode", "imlt_mode", "iary_num", "iotf_mode", "reserve1")
	head2 <- list(
		readChar(file_ptr, nchars=4),			# LC : P+O ロギングヘッダータイプ L0, L1, L2, ED
		readBin(file_ptr, integer(), size=4),	# LC : P+O レコード長
		readBin(file_ptr, numeric(), n=2, size=8),	# INF: P+O アンテナのオフセット値 (AZEL) [deg]
		readBin(file_ptr, numeric(), n=3, size=8),	# INF: P+O 副鏡のオフセット値 [mm] 0:X, 1:Z1, 2:Z2
		readBin(file_ptr, numeric(), n=6, size=8),	# SET: P+O 天体の座標値 [deg] [*][0]:RADEC, [*][1]:LB, [*][2]:AZEL
		readBin(file_ptr, numeric(), n=1, size=8),	# SET: P+O 観測者が与えた後退速度 [m/s]
		readBin(file_ptr, numeric(), n=1, size=8),	# SET: P+O マップのポジションアングル [deg]
		readBin(file_ptr, numeric(), n=1, size=8),	# SET: P+O (MLT) マルチビーム位置角オフセット値 [deg] (初期回転角と同意)
		readBin(file_ptr, numeric(), n=1, size=8),	# SET: --- スイッチング周波数 [Hz]
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# SET: --- 各アレイの周波数スイッチング間隔周波数 [Hz]
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# SET: P+O 各アレイの第一中間周波数 [Hz]
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# SET: P+O 各アレイの静止中心周波数 [Hz]
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# SET: P+O 各アレイの静止トラッキング周波数 [Hz] 
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# SET: P+O 各アレイのバックエンドのバンド幅 [Hz]; 2000^6, 1000^6, 500^6, 250^6, 125^6, 62.5^6, 31.25^6, 15.625^6
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# SET: P+O 周波数較正データ アレイの中心周波数 [Hz]
		readBin(file_ptr, numeric(), n=ARYMAX*2, size=8),	# SET: P+O 周波数較正データ 各測定点の周波数 [Hz]
		readBin(file_ptr, numeric(), n=ARYMAX*2, size=8),	# SET: P+O 周波数較正データ 各測定点のチャネル値
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# LC : P+O 各アレイのバックエンドのチャネル間隔 (観測バンド幅 / 出力ch数)
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# LC : P+O 各アレイのバックエンドの分解能 (観測バンド幅 / 出力ch数 * 1.00)
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# FL : P+O 各アレイの開口能率
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# FL : P+O 各アレイの主ビーム能率
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# FL : P+O 各アレイのアンテナ能率
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# FL : P+O 各アレイのFSS能率
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# FL : P+O 各アレイのアンテナ利得
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# FL : P+O 各アレイのビームサイズ
		readBin(file_ptr, numeric(), n=ARYMAX, size=8),	# FL : P+O 各アレイの直線偏波角度
		readBin(file_ptr, numeric(), n=BINMAX, size=8),	# SAM: P+O 1サイクル毎のビン毎の積分時間 [sec]
		readBin(file_ptr, numeric(), n=1, size=8),	# SET:   O (OTF) 1スキャンにおけるデータ取得間隔 [sec] 0:POSSW, 0.01-INTEG_TIME:OTF
		readBin(file_ptr, numeric(), n=15, size=8),	# ---: --- 予備領域
		readBin(file_ptr, integer(), n=ARYMAX, size=4),	# SET:   O (OTF) チャネルBIND値
		readBin(file_ptr, integer(), n=ARYMAX*2, size=4),	# SET:   O (OTF) 生データ保存チャネルレンジ
		readBin(file_ptr, integer(), n=1, size=4),		# SET: P+O キャリブ挿入のシーケンス間隔 [sec]
		readBin(file_ptr, integer(), n=ARYMAX, size=4),	# LC : P+O アレイのチャネル数 1-8192 (現状4096が最大)
		readBin(file_ptr, integer(), n=ARYMAX, size=4),	# FL : P+O 各アレイの強度較正体の温度
		readBin(file_ptr, integer(), n=ARYMAX, size=4),	# FL : P+O 各アレイの円偏波回転方向
		readBin(file_ptr, integer(), n=SUBARY* SPWMAX, size=4),	# SAM: P+O SPW-サブアレイ Channel averaging factor
		readBin(file_ptr, integer(), n=SUBARY* SPWMAX, size=4),	# SAM: P+O SPW-サブアレイ Start channel
		readBin(file_ptr, integer(), n=SUBARY* SPWMAX, size=4),	# SAM: P+O SPW-サブアレイ End channel
		readBin(file_ptr, numeric(), n=1, size=4),	# QL : P+O RMS自動停止下限値
		readBin(file_ptr, numeric(), n=10*2, size=4),	# QL : P+O POINTING ON点 AZEL値
		readBin(file_ptr, integer(), size=4),	# QL : P+O RMS自動判定モード
		readBin(file_ptr, integer(), size=4),	# QL : P+O RMS判定モード 
		readBin(file_ptr, integer(), size=4),	# QL : P+O RMS判定ON点番号
		readBin(file_ptr, integer(), n=65, size=4),	# ---: --- 予備領域
		readBin(file_ptr, raw(), n=1),	# SET: P+O 観測の属性 1:目的天体, 2:ポインティング
		readBin(file_ptr, raw(), n=1),	# SET: P+O 天体の座標系 0:RADEC, 1:LB, 2:AZEL
		readBin(file_ptr, raw(), n=1),	# SET: P+O スキャンオフセットの座標系 0:RADEC, 1:LB, 2:AZEL
		readBin(file_ptr, raw(), n=1),	# SET: P+O スイッチングモード 1:POS, 5:FREQ1, 6:FREQ2
		readBin(file_ptr, raw(), n=1),	# SET: P+O 後退速度の定義系 1:RADIO, 2:OPTICAL
		readBin(file_ptr, raw(), n=1),	# SET: P+O 後退速度の基準系 1:LSR, 2:HEL
		readBin(file_ptr, raw(), n=ARYMAX),	# SET: P+O (MLT) 各アレイのマルチビーム受信機チャネル番号 0:OFF, 1-16:マルチビーム観測
		readBin(file_ptr, raw(), n=ARYMAX),	# SET: P+O 各アレイの中間周波数減衰器の加算値 0-9 (Tsysと強度計算で使用)
		readBin(file_ptr, raw(), n=ARYMAX),	# SET: P+O 各アレイのリファレンスシンセサイザ番号 1-8
		readBin(file_ptr, raw(), n=ARYMAX),	# SET: P+O 各アレイの使用宣言フラグ 0:未使用, 1:使用
		readBin(file_ptr, raw(), n=1),	# SAM: P+O ビン数 1-BINMAX
		readBin(file_ptr, raw(), n=SUBARY),	# SAM: P+O SPW-サブアレイ Number of SPWectral windows 1-SPWMAX
		readChar(file_ptr, nchars=20), 		# MNG: P+O ロギングID YYYYMMDDHHMMSS
		readChar(file_ptr, nchars=16), 		# LC : P+O 観測開始年月日時分秒 (JST) YYYYMMDDHHMMSS
		readChar(file_ptr, nchars=20), 		# SET: P+O 観測指示書名
		readChar(file_ptr, nchars=16), 		# SET: P+O 観測天体名
		readChar(file_ptr, nchars=24), 		# SET: P+O 観測天体名のコメント
		readChar(file_ptr, nchars=40), 		# SET: P+O 観測者名
		readChar(file_ptr, nchars=10), 		# SET: P+O NEWSTARのグループ名
		readChar(file_ptr, nchars=10), 		# SET: P+O NEWSTARのプロジェクト名
		readChar(file_ptr, nchars=8), 		# SET: P+O 天体種別 RADEC, LB, AZEL, HOLO, SOLAR, COMET
		readChar(file_ptr, nchars=8), 		# SET: P+O 座標元期 B1950, J2000
		readChar(file_ptr, nchars=8), 		# SET: P+O スキャン種別 RADEC, LB, AZEL, SUNA, SUND
		readChar(file_ptr, nchars=120), 	# SET: P+O シーケンスパターン *:ON点, 1:OFF点
		readChar(file_ptr, nchars=ARYMAX*4), 	# SET: P+O 各アレイのサイドバンド種別 USB, LSB, DSB
		readChar(file_ptr, nchars=ARYMAX*10), 	# SET: P+O 各アレイの受信機フロントエンドの種別
		readChar(file_ptr, nchars=ARYMAX*2), 	# FL : P+O 各アレイのフィードホーン偏波の種類 R, L, H, V
		readChar(file_ptr, nchars=ARYMAX*5), 	# FL : P+O 各アレイの偏波の種類 CIRC, LINR
		readChar(file_ptr, nchars=4), 	# SAM: P+O サブアレイ分割設定 OFF, ON
		readChar(file_ptr, nchars=SUBARY* 6), 	# SAM: P+O 入力モード設定 XX, YY, XX&YY
		readChar(file_ptr, nchars=SUBARY* 64), 	# SAM: P+O BUTスケーリングファクタ設定 (ファイル名)
		readChar(file_ptr, nchars=SUBARY* 64), 	# SAM: P+O 再量子化再スケーリングファクタ設定 (ファイル名)
		readChar(file_ptr, nchars=SUBARY* 64), 	# SAM: P+O 再量子化 スケーリングファクタ設定 (ファイル名)
		readChar(file_ptr, nchars=SUBARY* 4), 	# SAM: P+O 高周波補正テーブル (逆重畳補正値) 作成 OFF, ON
		readChar(file_ptr, nchars=SUBARY* 4), 	# SAM: P+O 機能設定 HFS Cancellation OFF, ON
		readChar(file_ptr, nchars=SUBARY* 4), 	# SAM: P+O 機能設定 LFS Cancellation OFF, ON
		readChar(file_ptr, nchars=SUBARY* 7), 	# SAM: P+O 機能設定 Multiplexer selection BEFORE, AFTER
		readChar(file_ptr, nchars=SUBARY* SPWMAX* 4), 	# SAM: P+O SPW-サブアレイ Channel averaging mode OFF, ON
		readChar(file_ptr, nchars=SUBARY* SPWMAX* 16), 	# SAM: P+O SPW-サブアレイ Window function ID; NONE, HANNING, HAMMING, BARTLETT, BLACKMAN, BLACKMAN_HARRIS, WELCH
		readChar(file_ptr, nchars=4), 	# QL : P+O RMS判定アレイ名
		readBin(file_ptr, raw(), n=1),	# LC : P+O (MLT) 2beam観測のモード種別 0:2beam以外, 1:beam1, 2:beam2, 3:beam1&beam2 (TZ用)
		readBin(file_ptr, raw(), n=1),	# S&L: P+O (MLT) マルチビーム観測モード 0:OFF  2,4,6,8,16:ON (観測装置beam数, 4以上から新マルチ)
		readBin(file_ptr, raw(), n=1),	# LC : P+O (OTF) 観測アレイ数 1-ARYMAX
		readBin(file_ptr, raw(), n=1),	# SET: P+O (OTF) OTFモード 0:OFF, 1:ON
		readBin(file_ptr, raw(), n=509))	# ---: --- 予備領域
	names(head2) <- namelist
	return(head2)
}

#-------- Read scan records in SAM45 logging
readSAM45_data <- function( file_ptr ){
	namelist <- c("crec_type", "irec_len", "dant_prog", "dant_real", "dant_vel", "dmap_center", "dobs_pos", "dch_center", "dscan_pos_xy", "dscan_pos_xy_origin", "dsun_distance", "dsun_ps_angle", "dint_edmjd", "dint_stmjd", "dtsys", "dweather", "reserve8", "iinteg_time", "iline_no", "iscan_line", "iscan_no", "reserve4", "iary_no", "idmy_flag", "cary_name", "cint_edtm", "cint_sttm", "cscan_type",  "cscanid", "clog_time", "icdpc_endcode", "reserve1", "fary_data")
	log_rec <- list(
		readChar(file_ptr, nchars=4),			# LC : P+O ロギングヘッダータイプ L0, L1, L2, ED
		readBin(file_ptr, integer(), size=4),	# LC : P+O レコード長
		readBin(file_ptr, numeric(), n=2, size=8),	# INF: P   アンテナPROG値 (AZEL) [deg]
		readBin(file_ptr, numeric(), n=2, size=8),	# INF: P   アンテナREAL値 (AZEL) [deg]
		readBin(file_ptr, numeric(), size=8),	# INF: P   アンテナ視線速度 (VRAD) [m/s]
		readBin(file_ptr, numeric(), n=6, size=8),	# INF: P   マップセンター座標 [deg] [*][0]:RADEC, [*][1]:LB, [*][2]:AZEL
		readBin(file_ptr, numeric(), n=6, size=8),	# INF: P   観測点座標 [deg] [*][0]:RADEC, [*][1]:LB, [*][2]:AZEL
		readBin(file_ptr, numeric(), n=2, size=8),	# INF: P   マルチビームチャネルセンター座標 (RADEC) [deg]
		readBin(file_ptr, numeric(), n=2, size=8),	# INF: P   彗星 XY座標 [rad|km]
		readBin(file_ptr, numeric(), n=2, size=8),	# INF: P   彗星 スキャンオフセット [deg|km] 
		readBin(file_ptr, numeric(), size=8),	# INF: P   彗星 太陽距離 [km]
		readBin(file_ptr, numeric(), size=8),	# INF: P   彗星 太陽基準ポジションアングル [rad]
		readBin(file_ptr, numeric(), size=8),	# LC : P   積分終了時刻 (MJD)
		readBin(file_ptr, numeric(), size=8),	# LC : P   積分開始時刻 (MJD)
		readBin(file_ptr, numeric(), size=8),	# LC : P+O システム雑音温度 [K]
		readBin(file_ptr, numeric(), n=5, size=8),	# P+O 気象データ (取得は1[回/スキャン]); FL :     0:気温[℃], 1:気圧[hPa], 2:水蒸気圧[hPa], 3:風速[m/s], 4:風向(真北0度)
		readBin(file_ptr, numeric(), n=16, size=8),	# ---: --- 予備領域
		readBin(file_ptr, integer(), size=4),	# SET: P+O 積分時間 [sec]
		readBin(file_ptr, integer(), size=4),	# LC :   O (OTF) ON点のスキャン点カウンタ (初期値1)
		readBin(file_ptr, integer(), size=4),	# LC :   O (OTF) ON点のスキャン点カウンタ
		readBin(file_ptr, integer(), size=4),	# LC :   O (OTF) スキャントータル通番 0:ZERO, 1:R, 2:SKY, ...
		readBin(file_ptr, integer(), n=14, size=4),	# ---: --- 予備領域
		readBin(file_ptr, raw(), n=1),	# LC : P+O アレイの番号 41-72 (41 + アレイ番号(0-))
		readBin(file_ptr, raw(), n=1),	# LC : P+O スキャンデータ無効フラグ 0:有効, 1:無効
		readChar(file_ptr, nchars=4),	# LC : P+O アレイ名 A01-A32
		readChar(file_ptr, nchars=16),	# LC : P   積分終了時刻 (JST) YYYYMMDDHHMMSS
		readChar(file_ptr, nchars=32),	# LC : PSW 積分開始時刻 (JST) YYYYMMDDHHMMSS;LC : OTF 積分中心時刻 (JST) YYYYMMDDHHMMSS.sss
		readChar(file_ptr, nchars=6),	# LC : P+O スキャンタイプ PRE, ZERO, R, SKY, OFF, ON
		readChar(file_ptr, nchars=20),	# LC : P   スキャンID YYMMDDHHMMSS.99999
		readChar(file_ptr, nchars=16),	# LC : P   ロギング時刻 (JST) YYYYMMDDHHMMSS (移動により、2010/7以前はNULL)
		readChar(file_ptr, nchars=1),	# LC : P+O CDPCから受信した積分データに対する処理コード (調査用情報)
		readChar(file_ptr, nchars=111),	# ---: --- 予備領域
		readBin(file_ptr, numeric(), n=CHMAX, size=4))	# LC : P+O アレイ積分データ
	names(log_rec) <- namelist
	return(log_rec)
}

#-------- Integrated function to read SAM45 logging
readSAM45 <- function(fname){
	file_size <- file.info(fname)$size
	file_ptr <- file(fname, "rb")
	head1 <- readSAM45head1( file_ptr )
	head2 <- readSAM45head2( file_ptr )
	lineRec <- readSAM45_data(file_ptr)
	recNum <- (file_size - head1$irec_len - head2$irec_len - 24) / (lineRec$irec_len + 8)	# Calculate number of records
	SAM45spec <- lineRec$fary_data
	SAM45df <- data.frame(lineRec$cary_name, lineRec$cscan_type, lineRec$dint_stmjd, lineRec$dint_edmjd, 3600.0*(lineRec$dobs_pos[3]-lineRec$dmap_center[3])*cos(lineRec$dant_real[2]*DEGRAD), 3600.0*(lineRec$dobs_pos[6]-lineRec$dmap_center[6]), lineRec$dant_real[1], lineRec$dant_real[2], lineRec$dant_vel)
	for(rec_index in 2:recNum){
		lineRec <- readSAM45_data(file_ptr)
		SAM45spec <- rbind(SAM45spec, lineRec$fary_data)
		SAM45df <- rbind(SAM45df, data.frame(lineRec$cary_name, lineRec$cscan_type, lineRec$dint_stmjd, lineRec$dint_edmjd, 3600.0*(lineRec$dobs_pos[3]-lineRec$dmap_center[3])*cos(lineRec$dant_real[2]*DEGRAD), 3600.0*(lineRec$dobs_pos[6]-lineRec$dmap_center[6]), lineRec$dant_real[1], lineRec$dant_real[2], lineRec$dant_vel))
	}
	close(file_ptr)
	namelist <- c("cary_name", "cscan_type", "mjd_st", "mjd_ed", "dAZ", "dEL", "AZ", "EL", "Vrad")
	names(SAM45df) <- namelist
	return(list(head1, head2, SAM45spec, SAM45df))
}

#-------- Function to find Scan Pattern from SAM45 Log
scanPointing <- function(SAM45File){
	SAM45Log <- readSAM45(SAM45File)
	head1 <- SAM45Log[[1]]; head2 <- SAM45Log[[2]]; SAM45spec <- SAM45Log[[3]]; SAM45df <- SAM45Log[[4]]
	arrayIndex <- which(SAM45df$cary_name == SAM45df$cary_name[1])
	mjdRange <- min( SAM45df$mjd_st[arrayIndex] ):max( SAM45df$mjd_ed[arrayIndex] )
	#-------- Scan Type
	scanLen  <- length(mjdRange)
	ScanType <- AZ <- EL <- dAZ <- dEL <- Vrad <- rep( NA, scanLen )
	for(index in arrayIndex){
		timeIndex <- which( mjdRange >= SAM45df$mjd_st[index] & mjdRange <= SAM45df$mjd_ed[index] )
		ScanType[timeIndex] <- as.character(SAM45df$cscan_type[index])
		AZ[timeIndex] <- SAM45df$AZ[index]
		EL[timeIndex] <- SAM45df$EL[index]
		dAZ[timeIndex] <- SAM45df$dAZ[index]
		dEL[timeIndex] <- SAM45df$dEL[index]
		Vrad[timeIndex] <- SAM45df$Vrad[index]
	}
	scanDF <- data.frame(mjdRange, ScanType, AZ, EL, dAZ, dEL, Vrad)
	DF_label <- c('mjdSec', 'scanType', 'AZ', 'EL', 'dAZ', 'dEL', 'Vrad')
	names(scanDF) <- DF_label
	return(scanDF)
}
#-------- Function to find Scan Pattern from SAM45 Log
scanPattern <- function(SAM45File, prefix, IF_ID, threshFile){
    load(threshFile)
    SAM45Log <- readSAM45(SAM45File)
    head1 <- SAM45Log[[1]]; head2 <- SAM45Log[[2]]; SAM45spec <- SAM45Log[[3]]; SAM45df <- SAM45Log[[4]]
    arrayIndex <- which(SAM45df$cary_name == SAM45df$cary_name[1])
    mjdRange <- min( SAM45df$mjd_st[arrayIndex] ):max( SAM45df$mjd_ed[arrayIndex] )
    #-------- Scan Type
    scanLen  <- length(mjdRange)
    ScanType <- AZ <- EL <- dAZ <- dEL <- Vrad <- rep( NA, scanLen )
    for(index in arrayIndex){
        timeIndex <- which( mjdRange >= SAM45df$mjd_st[index] & mjdRange <= SAM45df$mjd_ed[index] )
        ScanType[timeIndex] <- as.character(SAM45df$cscan_type[index])
        AZ[timeIndex] <- SAM45df$AZ[index]
        EL[timeIndex] <- SAM45df$EL[index]
        dAZ[timeIndex] <- SAM45df$dAZ[index]
        dEL[timeIndex] <- SAM45df$dEL[index]
        Vrad[timeIndex] <- SAM45df$Vrad[index]
    }
	#-------- Covering PolariS Prefix
    prefix_index <- c()
    for(index in arrayIndex){
        prefix_index <- append(prefix_index, findPrefix(SAM45df$mjd_st[index], prefix))
        prefix_index <- append(prefix_index, findPrefix(SAM45df$mjd_ed[index], prefix))
    }
	prefix_index <- unique( prefix_index[prefix_index != -1] )
	PolarisFileNum <- length(prefix_index)
	#-------- MJD and Power in PolariS P file
	mjdSecPolaris <- numeric(0); powerIF <- list()
	for(IF_index in 1:length(IF_ID)){
		tempPower <- numeric(0)
		for(index in 1:PolarisFileNum){
			tempPower <- append( tempPower, bitThresh(sprintf("%s.P.%02d", prefix[prefix_index[index]], IF_ID[IF_index]), Thresh[,IF_index]))
			#tempPower <- append( tempPower, Apower(sprintf("%s.A.%02d", prefix[prefix_index[index]], IF_ID[IF_index])))
		}
		powerIF[[IF_index]] <- tempPower
	}
	mjdSecPolaris <- append(mjdSecPolaris, prefix2MJDsec(prefix[prefix_index[1]]) + seq(0, length(powerIF[[1]])-1, by=1))
	#-------- Matching between SAM45 and PolariS
	match_index <- rep(NA, scanLen)
	for( index in 1:scanLen){
        if( min(mjdSecPolaris) > mjdRange[index] - 1.5){ next }
        if( max(mjdSecPolaris) < mjdRange[index] - 0.5){ next }
        match_index[index] <- which(mjdSecPolaris > (mjdRange[index] - 1.5) & mjdSecPolaris <= (mjdRange[index] - 0.5))
    }
	#-------- Pack into a data frame
	scanDF <- data.frame(mjdRange, ScanType, AZ, EL, dAZ, dEL, Vrad); DF_label <- c('mjdSec', 'scanType', 'AZ', 'EL', 'dAZ', 'dEL', 'Vrad')
	for(IF_index in 1:length(IF_ID)){
		scanDF <- cbind( scanDF, powerIF[[IF_index]][match_index] )
		DF_label <- append(DF_label, sprintf('power%02d', IF_ID[IF_index]))
	}
	names(scanDF) <- DF_label
	return( list(scanDF=scanDF, head1=head1, head2=head2, SAM45spec=SAM45spec) )
}

#-------- 
bunchVec16 <- function(vec){ return(bunch_vec(vec, 16)) }

#-------- 
gauss4bit <- function(bitDist){ return(gaussNbit(bitDist, 16)) } 

#-------- Function to estimage power using 256-level histogram
bitPower <- function(fname){
	bitDist <- apply(readBitDist(fname), 2, bunchVec16)
	gaussResults <- apply(bitDist, 2, gauss4bit)
	return(1/gaussResults[1,]^2)
}

#-------- Function to estimate power using 256-level histogram
bitThresh <- function(fname, thresh){
	bitDist <- readBitDist(fname)
    nLevel <- length(thresh) + 1
    bunchLevel <- dim(bitDist)[1] / nLevel
    if(bunchLevel > 1){
        bunchVecNlevel <- function(vec){ return(bunch_vec(vec, bunchLevel)) }
        bitDist <- apply(bitDist, 2, bunchVecNlevel)
    }
	gThresh <- function(nsample){ return(gaussThresh(nsample, thresh)) }
	gaussResults <- apply(bitDist, 2, gThresh)
	return(1/gaussResults[1,]^2)
}
#-------- Function to estimate power using autocorrelation
Apower <- function(fname){
	A <- readPolariS(fname)
	chNum <- nrow(A); chRange <- floor(0.05*chNum):floor(0.95*chNum)
	return( apply(A[chRange,], 2, mean))
}
#-------- Function to calculate Tsys from Scan Pattern
scanTsys <- function(ScanDF, Tamb){
    nameList <- names(ScanDF)
    power_ptr <- grep('power', nameList); IFnum <- length(power_ptr)
    SAM45Files <- unique(ScanDF$FileName)
    newNameList <- nameList
    for(fileName in SAM45Files){
        tmpScan <- ScanDF[ScanDF$FileName == fileName,]
        tmpR_index <- which(tmpScan$scanType == 'R')
        #-------- Median Window Filter
        medR <- median( tmpScan$power00[tmpR_index]); sdR <- sd( tmpScan$power00[tmpR_index])
        R_index <- which(abs(tmpScan$power00 - medR) < 3.0*sdR)
        #-------- On and Off scans
        index <- which(tmpScan$scanType == 'OFF' | tmpScan$scanType == 'ON')
        OutOfR_index <- which( tmpScan$mjdSec > max(tmpScan$mjdSec[R_index]))
        for(IF_index in 1:IFnum){
            IF_ID <- as.integer(strsplit(nameList[power_ptr[IF_index]], "power")[[1]][2])
            Tsys   <- rep(NA, length(tmpScan$mjdSec))
            RPower <- predict(smooth.spline(tmpScan$mjdSec[R_index], tmpScan[[power_ptr[IF_index]]][R_index], spar=1.0), tmpScan$mjdSec)$y
            RPower[OutOfR_index] <- RPower[max(R_index)]
            Tsys[index]  <- Tamb / (RPower[index] / tmpScan[[power_ptr[IF_index]]][index] - 1.0)
            if(fileName == SAM45Files[1]){ newNameList <- append(newNameList, sprintf('Tsys%02d', IF_ID))}
            tmpScan <- cbind(tmpScan, Tsys)
        }
        names(tmpScan) <- newNameList
        if( fileName == SAM45Files[1]){
            resultDF <- tmpScan
        } else {
            resultDF <- rbind(resultDF, tmpScan)
        }
    }
    return(resultDF)
}

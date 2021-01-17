import crcmod
import copy
from selfdrive.car.hyundai.values import CAR, CHECKSUM

hyundai_checksum = crcmod.mkCrcFun(0x11D, initCrc=0xFD, rev=False, xorOut=0xdf)


def create_lkas11(packer, frame, car_fingerprint, apply_steer, steer_req,
                  lkas11, sys_warning, sys_state, CC, bus = 0 ):
 
  values = lkas11
  #values = copy.deepcopy( lkas11 )
  values["CF_Lkas_LdwsSysState"] = sys_state
  values["CF_Lkas_SysWarning"] = 3 if sys_warning else 0
  #values["CF_Lkas_LdwsLHWarning"] = left_lane_depart
  #values["CF_Lkas_LdwsRHWarning"] = right_lane_depart
  values["CR_Lkas_StrToqReq"] = apply_steer
  values["CF_Lkas_ActToi"] = steer_req
  values["CF_Lkas_MsgCount"] = frame % 0x10

  left_lane = CC.hudControl.leftLaneVisible
  right_lane = CC.hudControl.rightLaneVisible
  enabled = CC.enabled

  if car_fingerprint in [CAR.SONATA, CAR.PALISADE, CAR.KIA_NIRO_EV, CAR.SANTA_FE, CAR.IONIQ_EV_2020]:
    values["CF_Lkas_LdwsActivemode"] = int(left_lane) + (int(right_lane) << 1)
    values["CF_Lkas_LdwsOpt_USM"] = 2

    # FcwOpt_USM 5 = Orange blinking car + lanes
    # FcwOpt_USM 4 = Orange car + lanes
    # FcwOpt_USM 3 = Green blinking car + lanes
    # FcwOpt_USM 2 = Green car + lanes
    # FcwOpt_USM 1 = White car + lanes
    # FcwOpt_USM 0 = No car + lanes
    values["CF_Lkas_FcwOpt_USM"] = 2 if enabled else 1

    # SysWarning 4 = keep hands on wheel
    # SysWarning 5 = keep hands on wheel (red)
    # SysWarning 6 = keep hands on wheel (red) + beep
    # Note: the warning is hidden while the blinkers are on
    values["CF_Lkas_SysWarning"] = 4 if sys_warning else 0

  elif car_fingerprint == CAR.HYUNDAI_GENESIS:
    # This field is actually LdwsActivemode
    # Genesis and Optima fault when forwarding while engaged
    values["CF_Lkas_LdwsActivemode"] = 2
  elif car_fingerprint == CAR.KIA_OPTIMA:
    values["CF_Lkas_LdwsActivemode"] = 0

  dat = packer.make_can_msg("LKAS11", 0, values)[2]

  if car_fingerprint in CHECKSUM["crc8"]:
    # CRC Checksum as seen on 2019 Hyundai Santa Fe
    dat = dat[:6] + dat[7:8]
    checksum = hyundai_checksum(dat)
  elif car_fingerprint in CHECKSUM["6B"]:
    # Checksum of first 6 Bytes, as seen on 2018 Kia Sorento
    checksum = sum(dat[:6]) % 256
  else:
    # Checksum of first 6 Bytes and last Byte as seen on 2018 Kia Stinger
    checksum = (sum(dat[:6]) + dat[7]) % 256

  values["CF_Lkas_Chksum"] = checksum

  return packer.make_can_msg("LKAS11", bus, values)


def create_clu11(packer, frame, clu11, button, speed = None, bus = 0 ): 
  #values = copy.deepcopy( clu11 )
  values = clu11

  if speed != None:
    values["CF_Clu_Vanz"] = speed

  values["CF_Clu_CruiseSwState"] = button
  values["CF_Clu_AliveCnt1"] = frame % 0x10
  return packer.make_can_msg("CLU11", bus, values)


def create_lfa_mfa(packer, frame, enabled, Navi_SCC_Camera_Act, lfahda_mfc, hda_set_speed = 0):
  values = lfahda_mfc 
  if Navi_SCC_Camera_Act:  # 2:camera, 1:highway
    values["HDA_USM"] = 2
    values["HDA_SysWarning"] = 1 if enabled else 0
    values["HDA_Icon_State"] = 2 if enabled else 0   # 1:HDA(stanby),  2:HDA:white
    if Navi_SCC_Camera_Act == 2: # camera detect.
      values["HDA_Active"] = 1 if enabled else 0

  else:
    values["HDA_SysWarning"] = 2 if enabled else 0

    
  #values["LFA_Icon_State"]  = 2 if enabled else 0
  #values["HDA_SysWarning"] = 1 if enabled else 0
  #values["HDA_Active"] = 1 if enabled else 0
  #values["HDA_Icon_State"] = 2 if enabled else 0   # 1:HDA(stanby),  2:HDA:white

  if hda_set_speed:
    values["HDA_VSetReq"] = hda_set_speed

  """
  values = {
    #"LFA_USM": lfahda_mfc["LFA_USM"],
    #"LFA_USM": 2,
    "LFA_Icon_State": 2 if enabled else 0,
    #"LFA_SysWarning": 0,
    #"HDA_USM": lfahda_mfc["HDA_USM"],
    "HDA_USM": 2,
    #"HDA_Active": 1 if hda_set_speed else 0,
    "HDA_Active": 1 if enabled else 0,
    #"HDA_Icon_State": 2 if hda_set_speed else 0,
    "HDA_Icon_State": 2 if enabled else 0,
    "HDA_VSetReq": hda_set_speed,

    #"HDA_USM": 2,
    #"HDA_C_State": 5 if enabled else 0,
    #"HDA_VSetReq": hda_speed_limit if enabled else 0,
    #"LFA_SysWarning": 0,
    #"LFA_Icon_State": 2 if enabled else 0,
    #"LFA_USM": 2,
  }
  """
  # HDA_USM 2 = ?
  # HDA_VSetReq = HDA speed limit

  # LFA_SysWarning 1 = "Switching to HDA", short beep
  # LFA_SysWarning 2 = "Switching to Smart Cruise control", short beep
  # LFA_SysWarning 3 =  LFA error

  # LFA_Icon_State 0 = no wheel
  # LFA_Icon_State 1 = white wheel
  # LFA_Icon_State 2 = green wheel
  return packer.make_can_msg("LFAHDA_MFC", 0, values)



def create_mdps12(packer, frame, mdps12):
  #values = copy.deepcopy( mdps12 )
  values = mdps12
  values["CF_Mdps_ToiActive"] = 0
  values["CF_Mdps_ToiUnavail"] = 1
  values["CF_Mdps_MsgCount2"] = frame % 0x100
  values["CF_Mdps_Chksum2"] = 0

  dat = packer.make_can_msg("MDPS12", 2, values)[2]
  checksum = sum(dat) % 256
  values["CF_Mdps_Chksum2"] = checksum

  return packer.make_can_msg("MDPS12", 2, values)
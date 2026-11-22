# W16 reading program for the board (CircuitPython), v2. Copy this file to the CIRCUITPY drive as code.py.
#
# Needs in CIRCUITPY/lib: adafruit_tmag5273.py, adafruit_bus_device/, adafruit_register/
# (Adafruit_CircuitPython_TMAG5273 at commit 1e7e42b was read for every name used here.)
#
# Output over USB serial, one line every 0.1 s:
#   D,<ms since start>,<range mT: 266 or 133>,<x uT>,<y uT>,<z uT>,<temperature C>,<conversion status>,<device status>
# The range alternates: BLOCK lines at +-266 mT, then BLOCK lines at +-133 mT, with a pause after each switch.
# A settings line starting with "#" is printed at each switch.
# Temperature compensation is OFF (the die is not the magnets); the die temperature is on every line.
# NOT YET RUN ON HARDWARE.

import time

import board

import adafruit_tmag5273 as tm

BLOCK = 8          # data lines per range before switching
SETTLE = 0.4       # seconds to wait after a range switch before the first line

make_i2c = getattr(board, "STEMMA_I2C", None) or board.I2C
sensor = tm.TMAG5273(make_i2c())

# Set everything the measurement depends on, whatever the driver's defaults are.
sensor.magnetic_channels = tm.MAG_CH_XYZ
sensor.conversion_average = tm.CONV_AVG_32X
sensor.low_noise = True
sensor.mag_tempco = tm.TEMPCO_NONE
sensor.temperature_enabled = True
sensor.operating_mode = tm.MODE_CONTINUOUS


def set_range(wide):
    sensor.xy_range_wide = wide
    sensor.z_range_wide = wide
    time.sleep(SETTLE)
    print("# w16 code.py v2 id=0x%02X wide_xy=%d wide_z=%d avg=%d tempco=%d low_noise=%d mode=%d" % (
        sensor.device_id, sensor.xy_range_wide, sensor.z_range_wide, sensor.conversion_average,
        sensor.mag_tempco, sensor.low_noise, sensor.operating_mode))


t0 = time.monotonic()
wide = True
while True:
    set_range(wide)
    for _ in range(BLOCK):
        x, y, z = sensor.magnetic
        print("D,%d,%d,%.1f,%.1f,%.1f,%.2f,%d,%d" % (
            int((time.monotonic() - t0) * 1000), 266 if wide else 133, x, y, z, sensor.temperature,
            sensor.conversion_status, sensor.device_status))
        time.sleep(0.1)
    wide = not wide

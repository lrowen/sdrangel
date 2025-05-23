///////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2012 maintech GmbH, Otto-Hahn-Str. 15, 97204 Hoechberg, Germany //
// written by Christian Daniel                                                   //
// Copyright (C) 2015-2020, 2022 Edouard Griffiths, F4EXB <f4exb06@gmail.com>    //
// Copyright (C) 2020 Jon Beniston, M7RCE <jon@beniston.com>                     //
//                                                                               //
// This program is free software; you can redistribute it and/or modify          //
// it under the terms of the GNU General Public License as published by          //
// the Free Software Foundation as version 3 of the License, or                  //
// (at your option) any later version.                                           //
//                                                                               //
// This program is distributed in the hope that it will be useful,               //
// but WITHOUT ANY WARRANTY; without even the implied warranty of                //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                  //
// GNU General Public License V3 for more details.                               //
//                                                                               //
// You should have received a copy of the GNU General Public License             //
// along with this program. If not, see <http://www.gnu.org/licenses/>.          //
///////////////////////////////////////////////////////////////////////////////////

#ifndef PLUGINS_SAMPLESOURCE_USRPOUTPUT_USRPOUTPUTSETTINGS_H_
#define PLUGINS_SAMPLESOURCE_USRPOUTPUT_USRPOUTPUTSETTINGS_H_

#include <stdint.h>

#include <QByteArray>
#include <QString>

/**
 * These are the settings individual to each hardware channel or software Tx chain
 * Plus the settings to be saved in the presets
 */
struct USRPOutputSettings
{
    int      m_masterClockRate;
    // global settings to be saved
    uint64_t m_centerFrequency;
    int      m_devSampleRate;
    int      m_loOffset;
    // channel settings
    uint32_t m_log2SoftInterp;
    float    m_lpfBW;        //!< Analog lowpass filter bandwidth (Hz)
    uint32_t m_gain;         //!< Optimally distributed gain (dB)
    QString  m_antennaPath;
    QString  m_clockSource;
    bool     m_transverterMode;
    qint64   m_transverterDeltaFrequency;
    uint8_t  m_gpioDir;      //!< GPIO pin direction; 0 ATR (automatic transmit/receive), 1 output
    uint8_t  m_gpioPins;     //!< GPIO pins levels for outputs
    bool     m_useReverseAPI;
    QString  m_reverseAPIAddress;
    uint16_t m_reverseAPIPort;
    uint16_t m_reverseAPIDeviceIndex;

    USRPOutputSettings();
    void resetToDefaults();
    QByteArray serialize() const;
    bool deserialize(const QByteArray& data);
    void applySettings(const QStringList& settingsKeys, const USRPOutputSettings& settings);
    QString getDebugString(const QStringList& settingsKeys, bool force=false) const;
};

#endif /* PLUGINS_SAMPLESOURCE_USRPOUTPUT_USRPOUTPUTSETTINGS_H_ */

///////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2023 Jon Beniston, M7RCE <jon@beniston.com>                     //
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

#include <QDebug>

#include <complex.h>

#include "dsp/dspengine.h"
#include "dsp/fftfactory.h"
#include "util/db.h"

#include "freqscanner.h"
#include "freqscannersink.h"

FreqScannerSink::FreqScannerSink() :
        m_channel(nullptr),
        m_channelSampleRate(48000),
        m_channelFrequencyOffset(0),
        m_scannerSampleRate(33320),
        m_centerFrequency(0),
        m_messageQueueToChannel(nullptr),
        m_fftSequence(-1),
        m_fft(nullptr),
        m_fftCounter(0),
        m_fftSize(1024),
        m_binsPerChannel(16),
        m_averageCount(0)
{
    applySettings(m_settings, QStringList(), true);
    applyChannelSettings(m_channelSampleRate, m_channelFrequencyOffset, 16, 4, true);
}

FreqScannerSink::~FreqScannerSink()
{
    if (m_fftSequence >= 0)
    {
        FFTFactory* fftFactory = DSPEngine::instance()->getFFTFactory();
        fftFactory->releaseEngine(m_fftSize, false, m_fftSequence);
    }
}

void FreqScannerSink::feed(const SampleVector::const_iterator& begin, const SampleVector::const_iterator& end)
{
    Complex ci;

    for (SampleVector::const_iterator it = begin; it != end; ++it)
    {
        Complex c(it->real(), it->imag());
        c *= m_nco.nextIQ();

        if (m_interpolatorDistance < 1.0f) // interpolate
        {
            while (!m_interpolator.interpolate(&m_interpolatorDistanceRemain, c, &ci))
            {
                processOneSample(ci);
                m_interpolatorDistanceRemain += m_interpolatorDistance;
            }
        }
        else // decimate (and filter)
        {
            if (m_interpolator.decimate(&m_interpolatorDistanceRemain, c, &ci))
            {
                processOneSample(ci);
                m_interpolatorDistanceRemain += m_interpolatorDistance;
            }
        }
    }
}

void FreqScannerSink::processOneSample(Complex &ci)
{
    ci /= SDR_RX_SCALEF;

    m_fft->in()[m_fftCounter] = ci;
    m_fftCounter++;
    if (m_fftCounter == m_fftSize)
    {
        // Apply windowing function
        m_fftWindow.apply(m_fft->in());

        // Perform FFT
        m_fft->transform();

        // Reorder (so negative frequencies are first) and average
        int halfSize = m_fftSize / 2;
        for (int i = 0; i < halfSize; i++) {
            m_fftAverage.storeAndGetAvg(m_magSq[i], magSq(i + halfSize), i);
        }
        for (int i = 0; i < halfSize; i++) {
            m_fftAverage.storeAndGetAvg(m_magSq[i + halfSize], magSq(i), i + halfSize);
        }

        if (m_fftAverage.nextAverage())
        {
            // Send results to channel
            if (getMessageQueueToChannel() && (m_settings.m_channelBandwidth != 0) && (m_binsPerChannel != 0))
            {
                FreqScanner::MsgScanResult* msg = FreqScanner::MsgScanResult::create(m_fftStartTime);
                QList<FreqScanner::MsgScanResult::ScanResult>& results = msg->getScanResults();

                for (int i = 0; i < m_settings.m_frequencySettings.size(); i++)
                {
                    if (m_settings.m_frequencySettings[i].m_enabled)
                    {
                        qint64 frequency = m_settings.m_frequencySettings[i].m_frequency;
                        qint64 startFrequency = m_centerFrequency - m_scannerSampleRate / 2;
                        qint64 diff = frequency - startFrequency;
                        float binBW = m_scannerSampleRate / (float)m_fftSize;

                        // Ignore results in upper and lower 12.5%, as there may be aliasing here from half-band filters
                        if ((diff >= m_scannerSampleRate / 8) && (diff < m_scannerSampleRate * 7 / 8))
                        {
                            int bin = std::round(diff / binBW);
                            int channelBins;

                            if (m_settings.m_frequencySettings[i].m_channelBandwidth.isEmpty())
                            {
                                channelBins = m_binsPerChannel;
                            }
                            else
                            {
                                int channelBW = m_settings.getChannelBandwidth(&m_settings.m_frequencySettings[i]);
                                channelBins = m_fftSize / (m_scannerSampleRate / (float)channelBW);
                            }

                            // Calculate power at that frequency
                            Real power;
                            if (m_settings.m_measurement == FreqScannerSettings::PEAK) {
                                power = peakPower(bin, channelBins);
                            } else {
                                power = totalPower(bin, channelBins);
                            }

                            // Calculate voice activity level if using voice trigger
                            Real voiceLevel = 0.0;
                            if (m_settings.m_voiceSquelchType == FreqScannerSettings::VoiceLsb) {
                                voiceLevel = voiceActivityLevel(bin, channelBins, true);
                            } else if (m_settings.m_voiceSquelchType == FreqScannerSettings::VoiceUsb) {
                                voiceLevel = voiceActivityLevel(bin, channelBins, false);
                            }

                            //qDebug() << "startFrequency:" << startFrequency << "m_scannerSampleRate:" << m_scannerSampleRate << "m_centerFrequency:" << m_centerFrequency << "frequency" << frequency << "bin" << bin << "power" << power << "voiceLevel" << voiceLevel;
                            FreqScanner::MsgScanResult::ScanResult result = {frequency, power, voiceLevel};
                            results.append(result);
                        }
                    }
                }
                getMessageQueueToChannel()->push(msg);
            }
            m_averageCount = 0;
            m_fftStartTime = QDateTime::currentDateTime();
        }
        m_fftCounter = 0;
    }
}

// Calculate total power in a channel containing the specified bin (i.e. sums adjacent bins in the same channel)
Real FreqScannerSink::totalPower(int bin, int channelBins) const
{
    // Skip bin between halfway between channels
    // Then skip first and last bins, to avoid spectral leakage (particularly at DC)
    int startBin = bin - channelBins / 2 + 1 + 1;
    Real magSqSum = 0.0f;
    for (int i = 0; i < channelBins - 2 - 1; i++) {
        int idx = startBin + i;
        if ((idx < 0) || (idx >= m_fftSize)) {
            continue;
        }
        magSqSum += m_magSq[idx];
    }
    Real db = CalcDb::dbPower(magSqSum);
    return db;
}

// Calculate peak power in a channel containing the specified bin
Real FreqScannerSink::peakPower(int bin, int channelBins) const
{
    // Skip bin between halfway between channels
    // Then skip first and last bins, to avoid spectral leakage (particularly at DC)
    int startBin = bin - channelBins/2 + 1 + 1;
    Real maxMagSq = std::numeric_limits<Real>::min();
    for (int i = 0; i < channelBins - 2 - 1; i++)
    {
        int idx = startBin + i;
        if ((idx < 0) || (idx >= m_fftSize)) {
            continue;
        }
        //qDebug() << "idx:" << idx << "power:" << CalcDb::dbPower(m_magSq[idx]);
        maxMagSq = std::max(maxMagSq, m_magSq[idx]);
    }
    Real db = CalcDb::dbPower(maxMagSq);
    return db;
}

Real FreqScannerSink::magSq(int bin) const
{
    Complex c = m_fft->out()[bin];
    Real v = c.real() * c.real() + c.imag() * c.imag();
    Real magsq = v / (m_fftSize * m_fftSize);
    return magsq;
}

void FreqScannerSink::applyChannelSettings(int channelSampleRate, int channelFrequencyOffset, int scannerSampleRate, int fftSize, int binsPerChannel, bool force)
{
    qDebug() << "FreqScannerSink::applyChannelSettings:"
            << " channelSampleRate: " << channelSampleRate
            << " channelFrequencyOffset: " << channelFrequencyOffset
            << " scannerSampleRate: " << scannerSampleRate
            << " fftSize: " << fftSize
            << " binsPerChannel: " << binsPerChannel;

    if ((m_channelFrequencyOffset != channelFrequencyOffset) ||
        (m_channelSampleRate != channelSampleRate) || force)
    {
        m_nco.setFreq(-channelFrequencyOffset, channelSampleRate);
    }

    if ((m_channelSampleRate != channelSampleRate) || (m_scannerSampleRate != scannerSampleRate) || force)
    {
        m_interpolator.create(16, channelSampleRate, scannerSampleRate / 2.2); // Filter potential aliasing resulting from half-band filters
        m_interpolatorDistance = (Real) channelSampleRate / (Real)scannerSampleRate;
        m_interpolatorDistanceRemain = m_interpolatorDistance;
    }

    if ((m_fftSize != fftSize) || force)
    {
        FFTFactory* fftFactory = DSPEngine::instance()->getFFTFactory();
        if (m_fftSequence >= 0) {
            fftFactory->releaseEngine(fftSize, false, m_fftSequence);
        }
        m_fftSequence = fftFactory->getEngine(fftSize, false, &m_fft);
        m_fftCounter = 0;
        m_fftStartTime = QDateTime::currentDateTime();
        m_fftWindow.create(FFTWindow::Hanning, fftSize);

        int averages = m_settings.m_scanTime * scannerSampleRate / 2 / fftSize;
        m_fftAverage.resize(fftSize, averages);
        m_magSq.resize(fftSize);
    }

    m_channelSampleRate = channelSampleRate;
    m_channelFrequencyOffset = channelFrequencyOffset;
    m_scannerSampleRate = scannerSampleRate;
    m_fftSize = fftSize;
    m_binsPerChannel = binsPerChannel;
}

void FreqScannerSink::applySettings(const FreqScannerSettings& settings, const QStringList& settingsKeys, bool force)
{
    qDebug() << "FreqScannerSink::applySettings:"
             << settings.getDebugString(settingsKeys, force)
             << " force: " << force;

    if (settingsKeys.contains("scanTime") || force)
    {
        int averages = settings.m_scanTime * m_scannerSampleRate / 2 / m_fftSize;
        m_fftAverage.resize(m_fftSize, averages);
    }

    if (force) {
        m_settings = settings;
    } else {
        m_settings.applySettings(settingsKeys, settings);
    }
}

// Voice activity detection for SSB signals
// Detects voice by looking for formant-like structure (broad spectral peaks)
// Returns a value from 0.0 (no voice) to 1.0 (strong voice signature)
Real FreqScannerSink::voiceActivityLevel(int bin, int channelBins, bool isLSB) const
{
    // Voice band in SSB is typically 100-3000 Hz from carrier
    // We look for 2-4 formant peaks with bandwidth 50-200 Hz each

    int startBin = bin - channelBins / 2 + 1;
    int endBin = startBin + channelBins - 1;

    if (startBin < 0 || endBin >= m_fftSize) {
        return 0.0;
    }

    // Calculate bin bandwidth in Hz
    float binBW = m_scannerSampleRate / (float)m_fftSize;

    // For LSB, spectrum is reversed - flip the search direction
    int step = isLSB ? -1 : 1;
    int searchStart = isLSB ? endBin : startBin;
    int searchEnd = isLSB ? startBin : endBin;

    // Find peaks above noise floor
    QVector<int> peakBins;
    QVector<Real> peakMags;

    // Calculate average noise floor
    Real noiseFloor = 0.0;
    int noiseCount = 0;
    for (int i = startBin; i <= endBin; i++) {
        noiseFloor += m_magSq[i];
        noiseCount++;
    }
    noiseFloor = (noiseCount > 0) ? (noiseFloor / noiseCount) : 1e-12;
    Real threshold = noiseFloor * 3.0; // 4.77 dB above noise

    // Simple peak detection
    int i = searchStart;
    while ((isLSB && i >= searchEnd) || (!isLSB && i <= searchEnd))
    {
        if (m_magSq[i] > threshold)
        {
            // Found potential peak start
            Real peakMag = m_magSq[i];
            int peakBin = i;

            // Find local maximum
            i += step;
            while ((isLSB && i >= searchEnd) || (!isLSB && i <= searchEnd))
            {
                if (m_magSq[i] > peakMag) {
                    peakMag = m_magSq[i];
                    peakBin = i;
                    i += step;
                } else if (m_magSq[i] > threshold) {
                    i += step;
                } else {
                    break; // Peak ended
                }
            }

            peakBins.append(peakBin);
            peakMags.append(peakMag);
        }
        i += step;
    }

    if (peakBins.size() < 2) {
        return 0.0; // Need at least 2 peaks for voice
    }

    // Measure peak bandwidths
    int broadPeakCount = 0;
    for (int p = 0; p < peakBins.size(); p++)
    {
        int peakBin = peakBins[p];
        Real peakMag = peakMags[p];
        Real halfPower = peakMag * 0.5; // 3dB point

        // Measure bandwidth at half power (-3dB)
        int bwCount = 1; // Peak bin itself

        // Search left
        for (int j = peakBin - 1; j >= startBin; j--) {
            if (m_magSq[j] > halfPower) {
                bwCount++;
            } else {
                break;
            }
        }

        // Search right
        for (int j = peakBin + 1; j <= endBin; j++) {
            if (m_magSq[j] > halfPower) {
                bwCount++;
            } else {
                break;
            }
        }

        float bandwidth = bwCount * binBW;

        // Voice formants are typically 50-200 Hz wide
        // CW signals are <50 Hz wide
        if (bandwidth >= 50.0 && bandwidth <= 200.0) {
            broadPeakCount++;
        }
    }

    // Check formant spacing (voice formants are typically 500-1500 Hz apart)
    bool goodSpacing = false;
    if (broadPeakCount >= 2 && peakBins.size() >= 2)
    {
        for (int p = 0; p < peakBins.size() - 1; p++)
        {
            int spacing = std::abs(peakBins[p + 1] - peakBins[p]);
            float spacingHz = spacing * binBW;
            if (spacingHz >= 400.0 && spacingHz <= 1800.0) {
                goodSpacing = true;
                break;
            }
        }
    }

    // Calculate voice activity score
    // 2-4 broad peaks with good spacing = strong voice signature
    float score = 0.0;

    if (broadPeakCount >= 2)
    {
        // Base score from number of broad peaks
        score = std::min(broadPeakCount / 4.0f, 1.0f);

        // Boost if spacing is good
        if (goodSpacing) {
            score = std::min(score * 1.5f, 1.0f);
        }

        // Penalize if too many narrow peaks (likely CW or noise)
        int narrowPeakCount = peakBins.size() - broadPeakCount;
        if (narrowPeakCount > broadPeakCount) {
            score *= 0.5;
        }
    }

    return score;
}

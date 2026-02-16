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

        // Accumulate voice activity levels on individual FFT (before averaging)
        // This captures sharp formant structure better than averaged spectrum
        int freqCount = m_settings.m_frequencySettings.size();
        if (m_voiceLevelSum.size() != freqCount) {
            m_voiceLevelSum.resize(freqCount);
            m_voiceLevelCount.resize(freqCount);
            m_voiceLevelSum.fill(0.0);
            m_voiceLevelCount.fill(0);
        }

        for (int i = 0; i < freqCount; i++)
        {
            if (m_settings.m_frequencySettings[i].m_enabled)
            {
                qint64 frequency = m_settings.m_frequencySettings[i].m_frequency;
                qint64 startFrequency = m_centerFrequency - m_scannerSampleRate / 2;
                qint64 diff = frequency - startFrequency;
                float binBW = m_scannerSampleRate / (float)m_fftSize;

                if ((diff >= m_scannerSampleRate / 8) && (diff < m_scannerSampleRate * 7 / 8))
                {
                    int bin = std::round(diff / binBW);
                    int channelBins;
                    if (m_settings.m_frequencySettings[i].m_channelBandwidth.isEmpty()) {
                        channelBins = m_binsPerChannel;
                    } else {
                        int channelBW = m_settings.getChannelBandwidth(&m_settings.m_frequencySettings[i]);
                        channelBins = m_fftSize / (m_scannerSampleRate / (float)channelBW);
                    }

                    Real voiceLevel = 0.0;
                    if (m_settings.m_voiceSquelchType == FreqScannerSettings::VoiceLsb) {
                        voiceLevel = voiceActivityLevel(bin, channelBins, true);
                    } else if (m_settings.m_voiceSquelchType == FreqScannerSettings::VoiceUsb) {
                        voiceLevel = voiceActivityLevel(bin, channelBins, false);
                    }

                    if (voiceLevel > 0.0) {
                        m_voiceLevelSum[i] += voiceLevel;
                        m_voiceLevelCount[i]++;
                    }
                }
            }
        }

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
                            int bin = std::round(diff / binBW); // Bin corresponding to the frequency
                            int channelBins; // Number of bins in the channel containing the frequency.
                                             // This is either the default (m_binsPerChannel)
                                             // or calculated based on the channel bandwidth if specified in settings for this frequency

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

                            // Use averaged voice activity level from individual FFTs
                            Real voiceLevel = 0.0;
                            if ((m_settings.m_voiceSquelchType == FreqScannerSettings::VoiceLsb ||
                                 m_settings.m_voiceSquelchType == FreqScannerSettings::VoiceUsb) &&
                                m_voiceLevelCount[i] > 0) {
                                voiceLevel = m_voiceLevelSum[i] / m_voiceLevelCount[i];
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

            // Reset voice level accumulators for next averaging period
            m_voiceLevelSum.fill(0.0);
            m_voiceLevelCount.fill(0);
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

// Compute magSq from raw FFT output with reordering (negative frequencies first)
Real FreqScannerSink::magSqFromRawFFT(int bin) const
{
    // m_magSq is reordered: negative freqs first, then positive
    // m_fft->out() is in standard FFT order: DC, positive freqs, negative freqs
    int halfSize = m_fftSize / 2;
    int fftBin;
    if (bin < halfSize) {
        // Negative frequencies: map to second half of FFT output
        fftBin = bin + halfSize;
    } else {
        // Positive frequencies: map to first half of FFT output
        fftBin = bin - halfSize;
    }
    return magSq(fftBin);
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

        // Resize voice level accumulators to match frequency count
        int freqCount = m_settings.m_frequencySettings.size();
        m_voiceLevelSum.resize(freqCount);
        m_voiceLevelCount.resize(freqCount);
        m_voiceLevelSum.fill(0.0);
        m_voiceLevelCount.fill(0);
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

    // Calculate average noise floor from raw FFT
    Real noiseFloor = 0.0;
    int noiseCount = 0;
    for (int i = startBin; i <= endBin; i++) {
        noiseFloor += magSqFromRawFFT(i);
        noiseCount++;
    }
    noiseFloor = (noiseCount > 0) ? (noiseFloor / noiseCount) : 1e-12;
    Real threshold = noiseFloor * 3.0; // 4.77 dB above noise

    // Simple peak detection using raw FFT data
    int i = searchStart;
    while ((isLSB && i >= searchEnd) || (!isLSB && i <= searchEnd))
    {
        Real binMagSq = magSqFromRawFFT(i);
        if (binMagSq > threshold)
        {
            // Found potential peak start
            Real peakMag = binMagSq;
            int peakBin = i;

            // Find local maximum
            i += step;
            while ((isLSB && i >= searchEnd) || (!isLSB && i <= searchEnd))
            {
                Real nextMagSq = magSqFromRawFFT(i);
                if (nextMagSq > peakMag) {
                    peakMag = nextMagSq;
                    peakBin = i;
                    i += step;
                } else if (nextMagSq > threshold) {
                    i += step;
                } else {
                    break; // Peak ended
                }
            }

            // Calculate frequency offset from carrier
            // For USB: carrier is at startBin, voice extends upward
            // For LSB: carrier is at endBin, voice extends downward
            int carrierBin = isLSB ? endBin : startBin;
            float freqOffset = std::abs(peakBin - carrierBin) * binBW;

            // Only include peaks within SSB voice bandwidth (0-3000 Hz from carrier)
            if (freqOffset <= 3000.0) {
                peakBins.append(peakBin);
                peakMags.append(peakMag);
            }
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
            if (magSqFromRawFFT(j) > halfPower) {
                bwCount++;
            } else {
                break;
            }
        }

        // Search right
        for (int j = peakBin + 1; j <= endBin; j++) {
            if (magSqFromRawFFT(j) > halfPower) {
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
        // Detect fundamental frequency (f0) by looking for harmonic structure
        // Voice has harmonics at f0, 2*f0, 3*f0, etc. with formants modulating them
        // Typical f0: male 85-180 Hz, female 165-255 Hz
        int carrierBin = isLSB ? endBin : startBin;

        // Try different f0 candidates in typical voice range (80-300 Hz)
        int maxHarmonics = 0;

        for (int f0Hz = 80; f0Hz <= 300; f0Hz += 10)
        {
            int f0Bins = (int)(f0Hz / binBW);
            int harmonicCount = 0;

            // Check for harmonics up to 3000 Hz (SSB bandwidth limit)
            for (int h = 1; h <= 10; h++)
            {
                int harmonicBin = carrierBin + (isLSB ? -1 : 1) * (h * f0Bins);

                // Check if any detected peak is near this harmonic (within ±30 Hz tolerance)
                int tolerance = (int)(30.0 / binBW);
                for (int p = 0; p < peakBins.size(); p++)
                {
                    if (std::abs(peakBins[p] - harmonicBin) <= tolerance)
                    {
                        harmonicCount++;
                        break;
                    }
                }

                // Stop checking beyond 3 kHz
                if (h * f0Hz > 3000) {
                    break;
                }
            }

            if (harmonicCount > maxHarmonics)
            {
                maxHarmonics = harmonicCount;
            }
        }

        // Need at least 3 harmonics aligned to confirm voice pitch structure
        // If mistuned by 1 kHz, formants won't align with any harmonic series
        if (maxHarmonics < 3) {
            return 0.0;
        }

        // Base score from number of broad peaks
        score = std::min(broadPeakCount / 4.0f, 1.0f);

        // Boost if spacing is good
        if (goodSpacing) {
            score = std::min(score * 1.5f, 1.0f);
        }

        // Boost if strong harmonic structure (4+ harmonics)
        if (maxHarmonics >= 4) {
            score = std::min(score * 1.2f, 1.0f);
        }

        // Penalize if too many narrow peaks (likely CW or noise) 
        // => This condition is ALWAYS true as there are always many more peaks than broad peaks
        // int narrowPeakCount = peakBins.size() - broadPeakCount;
        // if (narrowPeakCount > broadPeakCount) {
        //     score *= 0.5;
        // }

        // Just strongly penalize if there are no broad peaks at all
        if (broadPeakCount == 0) {
            score *= 0.1;
        }
    }

    return score;
}

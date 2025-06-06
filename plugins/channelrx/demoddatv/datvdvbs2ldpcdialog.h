///////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2012 maintech GmbH, Otto-Hahn-Str. 15, 97204 Hoechberg, Germany //
// written by Christian Daniel                                                   //
// Copyright (C) 2015-2019, 2021 Edouard Griffiths, F4EXB <f4exb06@gmail.com>    //
// Copyright (C) 2015 John Greb <hexameron@spam.no>                              //
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

#ifndef SDRGUI_GUI_DATVDVBS2LDPCDIALOG_H_
#define SDRGUI_GUI_DATVDVBS2LDPCDIALOG_H_

#include <QDialog>

namespace Ui {
    class DatvDvbS2LdpcDialog;
}

class DatvDvbS2LdpcDialog : public QDialog {
    Q_OBJECT
public:
    explicit DatvDvbS2LdpcDialog(QWidget* parent = nullptr);
    ~DatvDvbS2LdpcDialog();

    void setMaxTrials(int maxTrials);
    int getMaxTrials() { return m_maxTrials; }

private:
    Ui::DatvDvbS2LdpcDialog* ui;
    int m_maxTrials;

private slots:
    void accept();
    void on_maxTrials_valueChanged(int value);
};

#endif /* SDRGUI_GUI_DATVDVBS2LDPCDIALOG_H_ */

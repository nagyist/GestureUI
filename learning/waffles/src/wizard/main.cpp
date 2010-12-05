// ----------------------------------------------------------------
// The contents of this file are distributed under the CC0 license.
// See http://creativecommons.org/publicdomain/zero/1.0/
// ----------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include "../GClasses/GApp.h"
#include "../GClasses/GMacros.h"
#include "../GClasses/GRand.h"
#include "../GClasses/GTime.h"
#include "../GSup/GWidgets.h"
#include "../GClasses/GThread.h"
#include "Gui.h"
#include "usage.h"
#include <exception>
#include <string>
#include <vector>
#include <iostream>

using namespace GClasses;
using std::cout;
using std::vector;
using std::string;

class WizardDialog;


class WizardController : public ControllerBase
{
protected:
	WizardDialog* m_pRootDialog;
	WizardDialog* m_pCurrentDialog;
	WizardDialog* m_pGoodbyeDialog;
	UsageNode* m_pRootNode;
	vector<UsageNode*> m_globals;
	UsageNode* m_pDefaultNode;
	bool m_autoNext;

public:
	WizardController();
	virtual ~WizardController();

	void RunModal();
	void doNext();
	void autoNext() { m_autoNext = true; }
	UsageNode* globalUsageNode(const char* name);
};



class WizardDialog : public GWidgetDialog
{
protected:
	WizardController* m_pController;
	WizardDialog* m_pParent;
	UsageNode* m_pNode;
	GWidgetTextButton* m_pNextButton;
	GWidgetGrid* m_pGrid;
	GWidgetGrid* m_pGrid2;
	vector<WizardDialog*> m_children;
	int m_selection;
	vector<UsageNode*> m_choices;
	vector<UsageNode*> m_selected;
	bool m_doAutoNext;

	enum part_type
	{
		pt_string,
		pt_struct,
		pt_loner,
	};

	enum mode
	{
		mode_choose_one,
		mode_options,
		mode_goodbye,
		mode_struct,
		mode_first_mave,
	};

	mode m_mode;

public:
	WizardDialog(WizardController* pController, WizardDialog* pParent, UsageNode* pNode, int w, int h)
	: GWidgetDialog(w, h, 0xff220022), m_pController(pController), m_pParent(pParent), m_pNode(pNode)
	{
		m_pGrid = NULL;
		m_pGrid2 = NULL;
		m_selection = 0;
		m_doAutoNext = false;
		if(!pNode)
		{
			m_mode = mode_goodbye;
			m_pNextButton = new GWidgetTextButton(this, 860, 620, 130, 48, "Exit");
			GWidgetTextLabel* pMsg = new GWidgetTextLabel(this, 100, 180, 800, 364, "Done! A command that will perform the specified task has been printed to stdout. You may now execute that command (by pasting it into a console window), or use it in a script (by pasting it into the script). Have a nice day!", 0xff8888aa, 0, 2.0f);
			pMsg->wrap();
			return;
		}
		m_pNextButton = new GWidgetTextButton(this, 860, 620, 130, 48, "Next");
		if(pNode->parts().size() == 2 && pNode->choices().size() == 1 && strcmp(pNode->parts()[1].c_str(), pNode->choices()[0]->tok()) == 0)
			pNode = pNode->choices()[0];
		if(pNode->choices().size() > 0 && ((pNode->parts().size() == 1 && pNode->tok()[0] == '[') || (pNode->parts().size() == 2 && pNode->tok()[0] != '[' && pNode->tok()[0] != '<' && pNode->parts()[1].c_str()[0] == '[')))
		{
			m_mode = mode_choose_one;
			new GWidgetTextLabel(this, 25, 25, 290, 24, "Please choose a value for", 0xff8888aa, 0, 2.0f);
			new GWidgetTextLabel(this, 315, 25, 300, 24, pNode->tok(), 0xffaa8888, 0, 2.0f);
			GWidgetTextLabel* pDescr = new GWidgetTextLabel(this, 25, 50, 950, 96, pNode->descr(), 0xffaa8888);
			pDescr->wrap();
			m_pGrid = new GWidgetGrid(this, 3, 20, 150, 970, 450, 0xff000000);
			m_pGrid->setHeaderHeight(0);
			m_pGrid->setRowHeight(80);
			m_pGrid->setColumnWidth(0, 24);
			m_pGrid->setColumnWidth(1, 250);
			m_pGrid->setColumnWidth(2, 676);
			GWidgetBulletHole* pFirstBullet = NULL;
			for(size_t i = 0; i < pNode->choices().size(); i++)
			{
				UsageNode* pChoice = pNode->choices()[i];
				GWidgetBulletHole* pHole = new GWidgetBulletHole(m_pGrid, 0, 0, 22, 22);
				if(i == 0)
					pFirstBullet = pHole;
				m_pGrid->setWidget(0, i, pHole);
				m_pGrid->setWidget(1, i, new GWidgetTextLabel(m_pGrid, 0, 0, 250, 24, pChoice->tok(), 0xff88aa88, 0, 2.0f));
				pDescr = new GWidgetTextLabel(m_pGrid, 0, 0, 670, 64, pChoice->descr(), 0xffaa8888);
				pDescr->wrap();
				m_pGrid->setWidget(2, i, pDescr);
			}
			if(pFirstBullet)
				pFirstBullet->setChecked(true);
		}
		else if(pNode->tok()[0] == '<')
		{
			m_mode = mode_options;
			new GWidgetTextLabel(this, 25, 25, 290, 24, "Specify optional stuff for", 0xff8888aa, 0, 2.0f);
			new GWidgetTextLabel(this, 315, 25, 200, 24, pNode->tok(), 0xffaa8888, 0, 2.0f);
			new GWidgetTextLabel(this, 20, 125, 200, 25, "Available options:", 0xff8888aa, 0, 2.0f);
			m_pGrid = new GWidgetGrid(this, 3, 20, 150, 710, 400, 0xff002200);
			m_pGrid->setHeaderHeight(0);
			m_pGrid->setRowHeight(80);
			m_pGrid->setColumnWidth(0, 25);
			m_pGrid->setColumnWidth(1, 200);
			m_pGrid->setColumnWidth(2, 469);
			new GWidgetTextLabel(this, 750, 125, 200, 25, "Selected options:", 0xff8888aa, 0, 2.0f);
			m_pGrid2 = new GWidgetGrid(this, 2, 750, 150, 241, 400, 0xff002200);
			m_pGrid2->setHeaderHeight(0);
			m_pGrid2->setRowHeight(64);
			m_pGrid2->setColumnWidth(0, 25);
			m_pGrid2->setColumnWidth(1, 200);
			for(size_t i = 0; i < pNode->choices().size(); i++)
				m_choices.push_back(pNode->choices()[i]);
			populateLists();
		}
		else if(pNode->parts().size() == 1 && pNode->tok()[0] != '[')
		{
			m_mode = mode_options;
			m_doAutoNext = true;
		}
		else if(countStructureParts(pNode) > 0)
		{
			m_mode = mode_struct;
			new GWidgetTextLabel(this, 25, 25, 290, 24, "Please provide values for...", 0xff8888aa);
			new GWidgetTextLabel(this, 315, 25, 300, 24, pNode->tok(), 0xffaa8888);
			GWidgetTextLabel* pDescr = new GWidgetTextLabel(this, 25, 50, 950, 96, pNode->descr(), 0xffaa8888);
			pDescr->wrap();
			m_pGrid = new GWidgetGrid(this, 3, 20, 150, 970, 450, 0xff000022);
			m_pGrid->setHeaderHeight(0);
			m_pGrid->setRowHeight(80);
			m_pGrid->setColumnWidth(0, 120);
			m_pGrid->setColumnWidth(1, 300);
			m_pGrid->setColumnWidth(2, 530);
			int row = 0;
			for(size_t i = 0; i < pNode->parts().size(); i++)
			{
				if(partType(pNode, i) == pt_struct)
				{
					string arg = pNode->parts()[i];
					UsageNode* pChoice = pNode->choice(arg.c_str());
					if(!pChoice)
						pChoice = globalUsageNode(arg.c_str());
					m_pGrid->setWidget(0, row, new GWidgetTextLabel(m_pGrid, 0, 0, 120, 24, pChoice->tok(), 0xff8888aa));
					GWidgetTextBox* pTB = new GWidgetTextBox(m_pGrid, 0, 0, 300, 24);
					m_pGrid->setWidget(1, row, pTB);
					if(row == 0)
						setFocusWidget(pTB);
					pDescr = new GWidgetTextLabel(m_pGrid, 0, 0, 530, 64, pChoice->descr(), 0xff8888aa);
					pDescr->wrap();
					m_pGrid->setWidget(2, row, pDescr);
					row++;
				}
			}
		}
		else
		{
			m_mode = mode_struct;
			m_doAutoNext = true;
		}
	}

	virtual ~WizardDialog()
	{
		clearChildren();
	}

	static bool isMultiUseFlag(const char* szFlag)
	{
		if(strcmp(szFlag, "-addlayer") == 0)
			return true;
		if(strcmp(szFlag, "[instance_count]") == 0)
			return true;
		return false;
	}

	void onShowDialog()
	{
		if(m_doAutoNext)
			m_pController->autoNext();
	}

	void populateLists()
	{
		m_pGrid->flushItems();
		m_pGrid2->flushItems();
		for(size_t i = 0; i < m_choices.size(); i++)
		{
			UsageNode* pChoice = m_choices[i];
			GWidgetBulletHole* pHole = new GWidgetBulletHole(m_pGrid, 0, 0, 22, 22);
			m_pGrid->setWidget(0, i, pHole);
			string s;
			pChoice->sig(&s);
			GWidgetTextLabel* pSig = new GWidgetTextLabel(m_pGrid, 0, 0, 200, 64, s.c_str(), 0xff00a0a0);
			pSig->wrap();
			m_pGrid->setWidget(1, i, pSig);
			GWidgetTextLabel* pDescr = new GWidgetTextLabel(m_pGrid, 0, 0, 469, 64, pChoice->descr(), 0xff00a0a0);
			pDescr->wrap();
			m_pGrid->setWidget(2, i, pDescr);
		}
		for(size_t i = 0; i < m_selected.size(); i++)
		{
			UsageNode* pChoice = m_selected[i];
			GWidgetBulletHole* pHole = new GWidgetBulletHole(m_pGrid2, 0, 0, 22, 22);
			m_pGrid2->setWidget(0, i, pHole);
			string s;
			pChoice->sig(&s);
			GWidgetTextLabel* pSig = new GWidgetTextLabel(m_pGrid2, 0, 0, 200, 64, s.c_str(), 0xff00a0a0);
			pSig->wrap();
			m_pGrid2->setWidget(1, i, pSig);
		}
		m_selection = 0;
	}

	part_type partType(UsageNode* pNode, int part)
	{
		const char* name = pNode->parts()[part].c_str();
		if(name[0] == '<')
			return pt_loner;
		if(name[0] != '[')
			return pt_string;
		UsageNode* pChoice = pNode->choice(name);
		if(!pChoice)
			pChoice = globalUsageNode(name);
		if(pChoice->choices().size() > 0)
			return pt_loner;
		else
			return pt_struct;
	}

	UsageNode* globalUsageNode(const char* name)
	{
		return m_pController->globalUsageNode(name);
	}

	int countStructureParts(UsageNode* pNode)
	{
		int count = 0;
		for(size_t i = 0; i < pNode->parts().size(); i++)
		{
			if(partType(pNode, i) == pt_struct)
				count++;
		}
		return count;
	}

	WizardDialog* parent()
	{
		return m_pParent;
	}

	vector<WizardDialog*>& children()
	{
		return m_children;
	}

	void clearChildren()
	{
		for(vector<WizardDialog*>::iterator it = m_children.begin(); it != m_children.end(); it++)
			delete(*it);
		m_children.clear();
	}

	virtual void onCheckBulletHole(GWidgetBulletHole* pBullet)
	{
		for(int i = 0; i < m_pGrid->rowCount(); i++)
		{
			GWidgetBulletHole* pBH = (GWidgetBulletHole*)m_pGrid->widget(0, i);
			if(pBH == pBullet)
				m_selection = i;
			else
				pBH->setChecked(false);
		}
		if(m_pGrid2)
		{
			for(int i = 0; i < m_pGrid2->rowCount(); i++)
			{
				GWidgetBulletHole* pBH = (GWidgetBulletHole*)m_pGrid2->widget(0, i);
				if(pBH == pBullet)
					m_selection = m_pGrid->rowCount() + i;
				else
					pBH->setChecked(false);
			}
		}

		if(m_mode == mode_options)
		{
			// Move the option to the other side
			if(m_selection < (int)m_choices.size())
			{
				UsageNode* pChoice = m_choices[m_selection];
				if(!isMultiUseFlag(m_choices[m_selection]->tok()))
					m_choices.erase(m_choices.begin() + m_selection);
				m_selected.push_back(pChoice);
				populateLists();
			}
			else
			{
				UsageNode* pChoice = m_selected[m_selection - m_choices.size()];
				m_selected.erase(m_selected.begin() + m_selection - m_choices.size());
				if(!isMultiUseFlag(m_selected[m_selection - m_choices.size()]->tok()))
					m_choices.push_back(pChoice);
				populateLists();
			}
		}
	}

	bool createChildDialogs()
	{
		clearChildren();
		if(m_mode == mode_choose_one)
		{
			if(m_selection < 0)
			{
				MessageBoxDialog dlg("You must select a choice first. (Click in the circle to the left of the choice you want.)");
				RunPopup(&dlg);
				return false;
			}
			UsageNode* pChoice = m_pNode->choices()[m_selection];
			m_children.push_back(new WizardDialog(m_pController, this, pChoice, m_rect.w, m_rect.h));
		}
		else if(m_mode == mode_struct)
		{
			// Check values
			for(int i = 0; m_pGrid && i < m_pGrid->rowCount(); i++)
			{
				GWidgetTextLabel* pLabel = (GWidgetTextLabel*)m_pGrid->widget(1, i);
				if(pLabel->text().size() < 1)
				{
					MessageBoxDialog dlg("One or more of the fields is still blank.");
					RunPopup(&dlg);
					return false;
				}
				bool quot = false;
				for(size_t j = 0; j < pLabel->text().size(); j++)
				{
					if(pLabel->text()[j] == '"')
						quot = !quot;
					else if(!quot && pLabel->text()[j] == ' ')
					{
						MessageBoxDialog dlg("Spaces are not allowed unless the field is wrapped in quotes.");
						RunPopup(&dlg);
						return false;
					}
				}
			}

			// Make the child dialogs
			for(size_t i = 0; i < m_pNode->parts().size(); i++)
			{
				if(partType(m_pNode, i) == pt_loner)
				{
					string arg = m_pNode->parts()[i];
					UsageNode* pChoice = m_pNode->choice(arg.c_str());
					if(!pChoice)
						pChoice = globalUsageNode(arg.c_str());
					m_children.push_back(new WizardDialog(m_pController, this, pChoice, m_rect.w, m_rect.h));
				}
			}
		}
		else if(m_mode == mode_options)
		{
			for(size_t i = 0; i < m_selected.size(); i++)
			{
				UsageNode* pChoice = m_selected[i];
				m_children.push_back(new WizardDialog(m_pController, this, pChoice, m_rect.w, m_rect.h));
			}
		}
		else
			ThrowError("Unrecognized mode");
		return true;
	}

	virtual void onReleaseTextButton(GWidgetTextButton* pButton)
	{
		if(pButton == m_pNextButton)
		{
			if(m_pNode)
				m_pController->doNext();
			else
				m_pController->quit();
		}
	}

	void printCommand()
	{
		if(m_mode == mode_choose_one)
		{
			GAssert(m_children.size() == 1); // expected one child dialog
			if(partType(m_pNode, 0) == pt_string)
				cout << m_pNode->tok() << " ";
			m_children[0]->printCommand();
		}
		else if(m_mode == mode_struct)
		{
			int stringPos = 0;
			int structPos = 0;
			int lonerPos = 0;
			for(size_t i = 0; i < m_pNode->parts().size(); i++)
			{
				part_type pt = partType(m_pNode, i);
				if(pt == pt_string)
				{
					cout << m_pNode->parts()[i] << " ";
					stringPos++;
				}
				else if(pt == pt_struct)
				{
					GWidgetTextLabel* pLabel = (GWidgetTextLabel*)m_pGrid->widget(1, structPos++);
					cout << pLabel->text() << " ";
				}
				else
				{
					GAssert(pt == pt_loner); // unexpected value
					m_children[lonerPos++]->printCommand();
				}
			}
			GAssert((size_t)(stringPos + structPos + lonerPos) == m_pNode->parts().size());
		}
		else if(m_mode == mode_options)
		{
			if(partType(m_pNode, 0) == pt_string)
				cout << m_pNode->tok() << " ";
			for(size_t i = 0; i < m_children.size(); i++)
				m_children[i]->printCommand();
		}
		else
			ThrowError("Unrecognized mode");
	}
};


// -------------------------------------------------------------------------------

class WizardView : public ViewBase
{
protected:
	WizardDialog* m_pDialog;

public:
	WizardView(WizardController* pController, WizardDialog* pInitialDialog)
	: ViewBase(), m_pDialog(pInitialDialog)
	{
	}

	virtual ~WizardView()
	{
	}

	void setDialog(WizardDialog* pDialog)
	{
		m_pDialog = pDialog;
		pDialog->onShowDialog();
		update();
	}

	virtual void onChar(char c)
	{
		m_pDialog->handleChar(c);
	}

	virtual void onMouseDown(int nButton, int x, int y)
	{
		m_pDialog->pressButton(nButton, x - m_screenRect.x, y - m_screenRect.y);
	}

	virtual void onMouseUp(int nButton, int x, int y)
	{
		m_pDialog->releaseButton(nButton);
	}

	virtual bool onMousePos(int x, int y)
	{
		return m_pDialog->handleMousePos(x - m_screenRect.x, y - m_screenRect.y);
	}

protected:
	virtual void draw(SDL_Surface *pScreen)
	{
		// Clear the screen
		SDL_FillRect(pScreen, NULL/*&r*/, 0x000000);

		// Draw the dialog
		blitImage(pScreen, m_screenRect.x, m_screenRect.y, m_pDialog->image());
	}
};


// -------------------------------------------------------------------------------


WizardController::WizardController()
: ControllerBase()
{
	m_pRootNode = makeMasterUsageTree();
	m_pRootDialog = new WizardDialog(this, NULL, m_pRootNode, 1010, 690);
	m_pGoodbyeDialog = NULL;
	m_pCurrentDialog = m_pRootDialog;
	m_pView = new WizardView(this, m_pCurrentDialog);
	m_pDefaultNode = new UsageNode("", "");
	m_globals.push_back(makeAlgorithmUsageTree());
	m_globals.push_back(makeNeighborUsageTree());
	m_autoNext = false;
}

WizardController::~WizardController()
{
	delete(m_pView);
	delete(m_pRootDialog);
	delete(m_pGoodbyeDialog);
	delete(m_pRootNode);
	delete(m_pDefaultNode);
	for(vector<UsageNode*>::iterator it = m_globals.begin(); it != m_globals.end(); it++)
		delete(*it);
}

UsageNode* WizardController::globalUsageNode(const char* name)
{
	for(vector<UsageNode*>::iterator it = m_globals.begin(); it != m_globals.end(); it++)
	{
		if(strcmp((*it)->tok(), name) == 0)
			return *it;
	}
	m_pDefaultNode->setTok(name);
	return m_pDefaultNode;
}

void WizardController::doNext()
{
	// If there are children, pick the first child
	if(!m_pCurrentDialog->createChildDialogs())
		return;
	WizardDialog* pNext = NULL;
	if(m_pCurrentDialog->children().size() > 0)
		pNext = m_pCurrentDialog->children()[0];

	if(!pNext)
	{
		WizardDialog* pPar = m_pCurrentDialog->parent();
		if(pPar)
		{
			// Find the current dialog
			size_t i;
			for(i = 0; i < pPar->children().size(); i++)
			{
				if(pPar->children()[i] == m_pCurrentDialog)
					break;
			}
			if(i >= pPar->children().size())
				ThrowError("internal error"); // failed to find current dialog

			// Pick the next sibling if there is one
			if(i + 1 < pPar->children().size())
				pNext = pPar->children()[i + 1];

			if(!pNext)
			{
				// Pick the next-sibling of the nearest ancestor with a next-sibling
				while(true)
				{
					WizardDialog* pParPar = pPar->parent();
					if(!pParPar)
						break;
					
					// Find the parent dialog
					size_t i;
					for(i = 0; i < pParPar->children().size(); i++)
					{
						if(pParPar->children()[i] == pPar)
							break;
					}
					if(i >= pParPar->children().size())
						ThrowError("internal error"); // failed to find pPar dialog

					// Pick the next sibling of pPar
					if(i + 1 < pParPar->children().size())
					{
						pNext = pParPar->children()[i + 1];
						break;
					}
					pPar = pParPar;
				}
			}
		}
	}
	if(pNext)
	{
		// Swap in the next dialog
		m_pCurrentDialog = pNext;
		((WizardView*)m_pView)->setDialog(m_pCurrentDialog);
	}
	else
	{
		// Print the command and terminate
		cout << "The command is:\n\n";
		m_pRootDialog->printCommand();
		cout << "\n";
		cout.flush();
		delete(m_pGoodbyeDialog);
		m_pGoodbyeDialog = new WizardDialog(this, NULL, NULL, 1010, 690);
		((WizardView*)m_pView)->setDialog(m_pGoodbyeDialog);
	}
}

void WizardController::RunModal()
{
	double timeOld = GTime::seconds();
	double time;
	double timeUpdate = 0;
	m_pView->update();
	while(m_bKeepRunning)
	{
		if(m_autoNext)
		{
			m_autoNext = false;
			doNext();
		}
		time = GTime::seconds();
		double delta = time - timeOld;
		if(handleEvents(delta) || time - timeUpdate > 5.0)
		{
			m_pView->update();
			timeUpdate = time;
		}
		else
			GThread::sleep(50);
		timeOld = time;
	}
}

int main(int argc, char *argv[])
{
	int nRet = 0;
	try
	{
		WizardController c;
		c.RunModal();
	}
	catch(const std::exception& e)
	{
		fprintf(stderr, "%s\n", e.what());
		nRet = 1;
	}

	return nRet;
}


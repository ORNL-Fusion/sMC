#include <Wt/WApplication>
#include <Wt/WBreak>
#include <Wt/WContainerWidget>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WText>
#include <Wt/WFileUpload>
#include <boost/filesystem.hpp>
#include <Wt/WTable>
#include <Wt/WInteractWidget>
#include <Wt/WSignalMapper>
#include <Wt/WIntValidator>
#include <Wt/WLineEdit>
#include <libconfig.h++>
#include <sstream>
#include <Wt/WDialog>
#include <Wt/WMessageBox>

#include "webFaceApp.h"

template <class T>
inline std::string to_string(const T &t) {
		std::stringstream ss;
		ss << t;
		return ss.str();
}

namespace bf = boost::filesystem;

webFaceApp::webFaceApp(const Wt::WEnvironment &env) : Wt::WApplication(env)
{
		this->useStyleSheet("smc.css");

		setTitle("sMC webFace by dlg0");

		logContainer_ = new Wt::WContainerWidget(root());
		cfgContainer_ = new Wt::WContainerWidget(root());

		logContainer_->setId("log");
		cfgContainer_->setId("cfg");

		// Existing Run 

		existingRunButton_ = new Wt::WPushButton("Choose existing run", root());
		existingRunButton_->clicked().connect(this, &webFaceApp::listRuns);

		uploadRunButton_ = new Wt::WPushButton("Upload run", root());
		uploadRunButton_->clicked().connect(this, &webFaceApp::uploadRun);

		newRunButton_ = new Wt::WPushButton("Create new run", root());
		newRunButton_->clicked().connect(this, &webFaceApp::newRun);

		// Setup cfg inputs and validators
		cfgContainer_->addWidget(new Wt::WText("nRow [16-1024]"));
		nRowLineEdit_ = new Wt::WLineEdit(cfgContainer_);
		Wt::WIntValidator *nRowValidator = new Wt::WIntValidator(16,1024);
		nRowLineEdit_->setValidator(nRowValidator);
		cfgContainer_->addWidget(new Wt::WBreak());

		cfgContainer_->addWidget(new Wt::WText("nCol [16-1024]"));
		nColLineEdit_ = new Wt::WLineEdit(cfgContainer_);
		Wt::WIntValidator *nColValidator = new Wt::WIntValidator(16,1024);
		nColLineEdit_->setValidator(nRowValidator);
		cfgContainer_->addWidget(new Wt::WBreak());

		cfgContainer_->addWidget(new Wt::WText("Atomic Z [1-118]"));
		ZLineEdit_ = new Wt::WLineEdit(cfgContainer_);
		Wt::WIntValidator *ZValidator = new Wt::WIntValidator(1,118,ZLineEdit_);
		cfgContainer_->addWidget(new Wt::WBreak());

		cfgContainer_->addWidget(new Wt::WText("Atomic amu [1-294]"));
		amuLineEdit_ = new Wt::WLineEdit(cfgContainer_);
		Wt::WIntValidator *amuValidator = new Wt::WIntValidator(1,294,amuLineEdit_);
		cfgContainer_->addWidget(new Wt::WBreak());

		cfgContainer_->addWidget(new Wt::WText("Eqdsk filename [relative path]"));
		eqdskLineEdit_ = new Wt::WLineEdit(cfgContainer_);
		eqdskValidator_ = new Wt::WValidator(eqdskLineEdit_);
		cfgContainer_->addWidget(new Wt::WBreak());

		root()->addWidget(new Wt::WBreak());
		
}

void webFaceApp::logEntry(const std::string &text)
{
		logEntry(text,"black");
}

void webFaceApp::logEntry(const std::string &text, const std::string &color)
{
		Wt::WText *tmpText = new Wt::WText(text);
		tmpText->decorationStyle().setForegroundColor(Wt::WColor(color));

		if(logContainer_->count()>0) {
			Wt::WWidget *topWidget = logContainer_->widget(0);
			logContainer_->insertBefore(new Wt::WBreak(),topWidget);
			topWidget = logContainer_->widget(0);
			logContainer_->insertBefore(tmpText,topWidget);
		}
		else
		{
			logContainer_->addWidget(tmpText);
		}
}

void webFaceApp::newRun()
{
		Wt::WDialog newRunDialog("Upload Run");
		new Wt::WText("Run name", newRunDialog.contents());
		new Wt::WLineEdit(newRunDialog.contents());
		new Wt::WBreak(newRunDialog.contents());

		Wt::WPushButton newRunButton("OK",newRunDialog.contents());
		newRunButton.clicked().connect(&newRunDialog, &Wt::WDialog::accept);
	
		newRunDialog.exec();
}

void webFaceApp::uploadRun()
{
		Wt::WDialog uploadRunDialog("Upload Run");
		new Wt::WText("Run name", uploadRunDialog.contents());
		new Wt::WLineEdit(uploadRunDialog.contents());
		new Wt::WBreak(uploadRunDialog.contents());

		uploadRunDialog.contents()->addWidget(new Wt::WText("Cfg file"));
		cfg_upload_ = new Wt::WFileUpload(uploadRunDialog.contents());
		cfg_upload_->setFileTextSize(40);
		Wt::WPushButton *cfg_uploadButton = new Wt::WPushButton("Upload",uploadRunDialog.contents());
		cfg_uploadButton->clicked().connect(cfg_upload_, &Wt::WFileUpload::upload);
		cfg_uploadButton->clicked().connect(cfg_uploadButton, &Wt::WPushButton::disable);
		spoolFileLocation_ = new Wt::WText(uploadRunDialog.contents());
		cfg_upload_->uploaded().connect(this, &webFaceApp::fileUploaded);
		new Wt::WBreak(uploadRunDialog.contents());

		uploadRunDialog.contents()->addWidget(new Wt::WText("Eqdsk file"));
		Wt::WFileUpload *eqdsk_upload_ = new Wt::WFileUpload(uploadRunDialog.contents());
		Wt::WPushButton *eqdsk_uploadButton = new Wt::WPushButton("Upload",uploadRunDialog.contents());
		new Wt::WBreak(uploadRunDialog.contents());

		Wt::WPushButton uploadReadyButton("OK",uploadRunDialog.contents());
		uploadReadyButton.clicked().connect(&uploadRunDialog, &Wt::WDialog::accept);
	
		uploadRunDialog.exec();
}

void webFaceApp::listRuns()
{

		// Build the run selection dialog
	
		Wt::WDialog existingRunDialog("Choose Run");
		//runContainer_ = new Wt::WContainerWidget(existingRunDialog_->contents());
		//runContainer_->setId("runList");

		// List runs in scratch
		root()->addWidget(new Wt::WBreak());
		bf::path scratch ("/home/dg6/scratch/sMC");
		Wt::WSignalMapper<Wt::WText *> *runMap = new Wt::WSignalMapper<Wt::WText *>(this);
		runMap->mapped().connect(this, &webFaceApp::selectRun);
		std::vector<Wt::WText *> runs;

		if(bf::exists(scratch)) {
				if(bf::is_directory(scratch)) {

					logEntry("Scratch directory is good.");

					typedef std::vector<bf::path> vec;
					vec v;
					copy(bf::directory_iterator(scratch),
									bf::directory_iterator(), back_inserter(v));
					for(int i=0;i<v.size();i++) {
						runs.push_back(new Wt::WText(v[i].string()));
						existingRunDialog.contents()->addWidget(runs[i]);
						existingRunDialog.contents()->addWidget(new Wt::WBreak());
						//runContainer_->addWidget(runs[i]);
						//runContainer_->addWidget(new Wt::WBreak());
						runMap->mapConnect(runs[i]->clicked(),runs[i]);
					}
				}
		}
		else
		{
				logEntry("Scratch directory does not exist.");
		}

		Wt::WPushButton runSelectedButton("Run Selected",existingRunDialog.contents());
		runSelectedButton.clicked().connect(&existingRunDialog, &Wt::WDialog::accept);
		existingRunDialog.exec();

}


void webFaceApp::selectRun(Wt::WText *text)
{
		// Add output
		logEntry(text->text().narrow());

		// Search for .cfg file
		std::string cfgName = text->text().narrow()+"/sMC.cfg";
		bf::path cfgPath (cfgName);
		bool cfgFound = bf::exists(cfgPath);
		if(cfgFound) {

			logEntry("Found sMC.cfg file","blue");

			// Readin cfg file
			libconfig::Config cfg;
			cfg.readFile(cfgName.c_str());
			int _Z = cfg.lookup("Z");
			int amu = cfg.lookup("amu");
			unsigned int nRow = cfg.lookup("nRow");
			unsigned int nCol = cfg.lookup("nCol");
			std::string eqdskName = cfg.lookup("eqdskPath");
			bf::path eqdskPath(text->text().narrow());
			eqdskPath/=eqdskName;

			// Update form
			nRowLineEdit_->setText(to_string(nRow));
			nColLineEdit_->setText(to_string(nCol));
			ZLineEdit_->setText(to_string(_Z));
			amuLineEdit_->setText(to_string(amu));
			eqdskLineEdit_->setText(eqdskName);
			if(!bf::exists(eqdskPath)) eqdskLineEdit_->setStyleClass("Wt-invalid");
			std::cout << eqdskPath.string();
		}
		else
		{
			logEntry("ERROR: Could not find sMC.cfg file","red");
		}
}

void webFaceApp::fileUploaded()
{
	spoolFileLocation_->setText(cfg_upload_->spoolFileName());
	cfg_upload_->stealSpooledFile();
}

Wt::WApplication *createApplication(const Wt::WEnvironment &env)
{
		return new webFaceApp(env);
}

int main (int argc, char **argv)
{
		return WRun(argc, argv, &createApplication);
}

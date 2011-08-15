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

namespace bf = boost::filesystem;

class HelloApp : public Wt::WApplication
{
		public:
			HelloApp(const Wt::WEnvironment &env);
		private:
			Wt::WLineEdit *nameEdit_;
			Wt::WText *greeting_, *spoolFileLocation_;
			Wt::WFileUpload *cfg_upload_;
			Wt::WContainerWidget *logContainer_, *runContainer_;

			void greet();
			void fileUploaded();
			void selectRun(Wt::WText *text);
			void logEntry(const std::string &text);
			void logEntryColored(const std::string &text, const std::string &color);
};

HelloApp::HelloApp(const Wt::WEnvironment &env) : Wt::WApplication(env)
{
		this->useStyleSheet("smc.css");

		setTitle("sMC webFace by dlg0");

		logContainer_ = new Wt::WContainerWidget(root());
		runContainer_ = new Wt::WContainerWidget(root());

		logContainer_->setId("log");

		root()->addWidget(new Wt::WBreak());
		
		//greeting_ = new Wt::WText(root());

		//b->clicked().connect(this, &HelloApp::greet);
		//nameEdit_->enterPressed().connect(this, &HelloApp::greet);

		// Cfg file upload
		root()->addWidget(new Wt::WText("Upload cfg file"));
		cfg_upload_ = new Wt::WFileUpload(root());
		cfg_upload_->setFileTextSize(40);
		Wt::WPushButton *cfg_uploadButton = new Wt::WPushButton("Upload",root());
		cfg_uploadButton->clicked().connect(cfg_upload_, &Wt::WFileUpload::upload);
		cfg_uploadButton->clicked().connect(cfg_uploadButton, &Wt::WPushButton::disable);
		spoolFileLocation_ = new Wt::WText(root());
		cfg_upload_->uploaded().connect(this, &HelloApp::fileUploaded);

		// List runs in scratch
		root()->addWidget(new Wt::WBreak());
		bf::path scratch ("/home/dg6/scratch/sMC");
		Wt::WTable *runTable = new Wt::WTable(runContainer_);
		Wt::WSignalMapper<Wt::WText *> *runMap = new Wt::WSignalMapper<Wt::WText *>(this);
		runMap->mapped().connect(this, &HelloApp::selectRun);
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
						runContainer_->addWidget(runs[i]);
						runContainer_->addWidget(new Wt::WBreak());
						runMap->mapConnect(runs[i]->clicked(),runs[i]);
					}
				}
		}
		else
		{
				logEntry("Scratch directory does not exist.");
		}

}

void HelloApp::logEntry(const std::string &text)
{
		logEntryColored(text,"black");
}

void HelloApp::logEntryColored(const std::string &text, const std::string &color)
{
		Wt::WText *tmpText = new Wt::WText(text);
		tmpText->decorationStyle().setForegroundColor(Wt::WColor(color));
		logContainer_->addWidget(tmpText);
		logContainer_->addWidget(new Wt::WBreak());
}

void HelloApp::selectRun(Wt::WText *text)
{
		// Add output
		logEntry(text->text().narrow());

		// Search for .cfg file
		bf::path cfgPath (text->text().narrow()+"/sMC.cfg");
		bool cfgFound = bf::exists(cfgPath);
		if(cfgFound) {
			logEntry("Found sMC.cfg file :)");
		}
		else
		{
			logEntryColored("ERROR: Could not find sMC.cfg file","red");
		}
}

void HelloApp::fileUploaded()
{
	spoolFileLocation_->setText(cfg_upload_->spoolFileName());
	cfg_upload_->stealSpooledFile();
}

void HelloApp::greet()
{
		greeting_->setText("Hello "+ nameEdit_->text());
}

Wt::WApplication *createApplication(const Wt::WEnvironment &env)
{
		return new HelloApp(env);
}

int main (int argc, char **argv)
{
		return WRun(argc, argv, &createApplication);
}

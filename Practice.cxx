#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkActor.h>
#include <vtkNamedColors.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkAppendFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkAssembly.h>
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkTransform.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <random>

namespace Atom {
	/**
	 * Create Particle.
	 *
	 * @param center - The center of particle.
	 * @param radius - The radius of particle.
	 * @param color - The color of particle.
	 *
	 * @return The nucleas.
	 */
	vtkSmartPointer<vtkActor> CreateParticle(double* center, double radius, double* color) {
		vtkSmartPointer<vtkSphereSource> particle = vtkSmartPointer<vtkSphereSource>::New();
		particle->SetRadius(radius);
		particle->SetPhiResolution(100);
		particle->SetThetaResolution(100);
		particle->SetCenter(center[0], center[1], center[2]);

		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(particle->GetOutputPort());

		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(color[0], color[1], color[2]);
		return actor;
	}

	/**
	 * Create Nucleas with protons and neutrons inside it.
	 *
	 * @param renderer - Renderer.
	 * @param nucleasRadius - The radius of nucleas.
	 * @param protons - The number of protons.
	 * @param neutrons - The number of neutrons.
	 *
	 * @return The nucleas.
	 */
	vtkSmartPointer<vtkAssembly> CreateNucleas(vtkSmartPointer<vtkAssembly> nucleusAssembly, int nucleasRadius, int protons, int neutrons, vtkSmartPointer<vtkNamedColors> colors) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> dist(0, 1);
		std::uniform_real_distribution<double> theta_dist(0.0, 2.0 * M_PI);
		std::uniform_real_distribution<double> phi_dist(0.0, M_PI);

		vtkSmartPointer<vtkAppendFilter> nucleas = vtkSmartPointer<vtkAppendFilter>::New();

		int totalParticles = protons + neutrons;
		for (int i = 1; i <= totalParticles; i++) {
			// Calculate random coordinates with respect to virtual sphere radius
			double u = dist(gen);
			double R = nucleasRadius * cbrt(u);  // cbrt(u) to ensure uniform distribution within the volume

			double theta = theta_dist(gen);
			double phi = phi_dist(gen);
			double center[3] { R * sin(phi) * cos(theta),
							   R * sin(phi) * sin(theta), 
				               R * cos(phi) };
			
			double *color = (i <= protons) ? colors->GetColor3d("Proton").GetData() : colors->GetColor3d("Neutron").GetData();
			vtkSmartPointer<vtkActor> actor = CreateParticle(center, 0.9, color);
			nucleusAssembly->AddPart(actor);
		}

		return nucleusAssembly;
	}

	void RotateNucleas(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData) {
		vtkRenderWindowInteractor *interactor = static_cast<vtkRenderWindowInteractor*>(caller);
		vtkAssembly *assembly = static_cast<vtkAssembly*>(clientData);

		vtkTransform* transform = vtkTransform::SafeDownCast(assembly->GetUserTransform());
		if (!transform) {
			transform = vtkTransform::New();
			assembly->SetUserTransform(transform);
		}

		transform->RotateY(1.0);
		interactor->GetRenderWindow()->Render();
	}
}

int main(int, char* [])
{
	vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
	colors->SetColor("Proton", "Yellow");
	colors->SetColor("Neutron", "Blue");
	colors->SetColor("Electron", "Black");

	vtkSmartPointer<vtkAssembly> nucleasAssembly = vtkSmartPointer<vtkAssembly>::New();
	nucleasAssembly = Atom::CreateNucleas(nucleasAssembly, 5, 94, 144, colors);
	
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(nucleasAssembly);
	renderer->SetBackground(colors->GetColor3d("White").GetData());

	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetWindowName("Plutonium-238 Isotope");
	renderWindow->AddRenderer(renderer);

	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkNew<vtkInteractorStyleSwitch> styleSwitch;
	styleSwitch->SetCurrentStyleToTrackballCamera();
	renderWindowInteractor->SetInteractorStyle(styleSwitch);

	// Orientation Marker
	vtkNew<vtkOrientationMarkerWidget> markerWidget;
	double rgba[4]{ 0.0, 0.0, 0.0, 0.0 };
	colors->GetColor("Carrot", rgba);
	markerWidget->SetOutlineColor(rgba[0], rgba[1], rgba[2]);

	vtkNew<vtkAxesActor> marker;
	markerWidget->SetOrientationMarker(marker);

	markerWidget->SetInteractor(renderWindowInteractor);
	markerWidget->SetViewport(0.85, 0.75, 1, 1);
	markerWidget->EnabledOn();
	markerWidget->InteractiveOn();

	vtkCamera* camera = renderer->GetActiveCamera();
	camera->SetPosition(0.0, 0.0, 100.0);

	vtkSmartPointer<vtkCallbackCommand> nucleasCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	nucleasCallback->SetCallback(Atom::RotateNucleas);
	nucleasCallback->SetClientData(nucleasAssembly);
	renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, nucleasCallback);

	//renderWindow->SetFullScreen(true);
	renderWindow->SetSize(800, 600);
	renderWindow->Render();
	renderWindowInteractor->Initialize();
	renderWindowInteractor->CreateRepeatingTimer(1);
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}
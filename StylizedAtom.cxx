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
#include <vtkProperty2D.h>
#include <vtkAssembly.h>
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkTransform.h>
#include <vtkEllipseArcSource.h>
#include <vtkRegularPolygonSource.h>
#include <vtkMath.h>
#include <vtkPropCollection.h>
#include <vtkProperty.h>
#include <vtkLight.h>
#include <vtkLegendBoxActor.h>
#include <vtkTextProperty.h>
#include <vtkVersion.h>
#include <random>

#if VTK_VERSION_NUMBER >= 90020210809ULL
#define VTK_HAS_COW 1
#endif

#if VTK_HAS_COW
#include <vtkCameraOrientationWidget.h>
#endif

namespace Atom {

	/**
	 * Create a single particle source.
	 *
	 * @param center - The center of particle
	 * @param radius - The radius of particle
	 * @param color - The color of particle
	 *
	 * @return particleSource
	 */
	vtkSmartPointer<vtkSphereSource> CreateParticleSource(double* center, double radius, double* color);

	/**
	 * Create a single particle actor.
	 *
	 * @param center - The center of particle
	 * @param radius - The radius of particle
	 * @param color - The color of particle
	 *
	 * @return particleActor
	 */
	vtkSmartPointer<vtkActor> CreateParticle(double* center, double radius, double* color);

	/**
	 * Create a nucleas with protons and neutrons inside it.
	 *
	 * @param nucleusAssembly - Group of protons and neutrons as a single assembly.
	 * @param nucleasRadius - The radius of resulting nucleas
	 * @param protons - The number of protons
	 * @param neutrons - The number of neutrons
	 * @param color - The color of particle inside nucleas
	 * 
	 * @return nucleusAssembly
	 */
	vtkSmartPointer<vtkAssembly> CreateNucleas(vtkSmartPointer<vtkAssembly> nucleusAssembly,
                                               int nucleasRadius,
                                               int protons,
                                               int neutrons,
                                               vtkSmartPointer<vtkNamedColors> colors);

	/**
	 * Create electron shells (energy levels) around the nucleas.
	 *
	 * @param renderer - The renderer
	 * @param shellCollection - Collection of shell assemblies
	 * @param shells - The number of shells
	 * @param color - The color of particle inside nucleas
	 *
	 * @return shellCollection
	 */
	vtkSmartPointer<vtkPropCollection> CreateElectronShell(vtkSmartPointer<vtkRenderer> renderer,
                                                           vtkSmartPointer<vtkPropCollection> shellCollection,
                                                           int shells,
                                                           vtkSmartPointer<vtkNamedColors> colors);

	/**
	 * Apply rotation transformation to the nucleas.
	 *
	 * @param caller - Caller invoking the event
	 * @param eventId - event id
	 * @param clientData - The data to apply rotation
	 * @param callData - Callback data sent with the event invokation
	 *
	 * @return void
	 */
	void RotateNucleas(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

	/**
	 * Apply rotation transformation to each individual shell.
	 *
	 * @param caller - Caller invoking the event
	 * @param eventId - event id
	 * @param clientData - The data to apply rotation
	 * @param callData - Callback data sent with the event invokation
	 *
	 * @return void
	 */
	void RotateShell(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

	/**
	 * Adjust camera zooming
	 *
	 * @param caller - Caller invoking the event
	 * @param eventId - event id
	 * @param clientData - The data to apply rotation
	 * @param callData - Callback data sent with the event invokation
	 *
	 * @return void
	 */
	void AdjustZoom(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

	/**
	 * Create legend
	 *
	 * @param color - The color for legend text properties
	 * 
	 * @return vtkLegendBoxActor
	 */
	vtkSmartPointer<vtkLegendBoxActor> CreateLegendWidget(vtkSmartPointer<vtkNamedColors> colors);

	unsigned long dollyObserverTag;
	double currentZoomInFactor;
	const double ZOOM_IN_DELTA = 0.02;
} // namespace

int main(int, char* [])
{
	// Define custom colors
	vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
	colors->SetColor("Proton", "Red");
	colors->SetColor("Neutron", "Blue");
	colors->SetColor("Electron", "Coral");
	colors->SetColor("Shell", "Coral");

	// Create a nucleus assembly consisting of protons and neutrons.
	vtkSmartPointer<vtkAssembly> nucleasAssembly = vtkSmartPointer<vtkAssembly>::New();
	nucleasAssembly = Atom::CreateNucleas(nucleasAssembly, 4, 50, 50, colors);

	// Add Nucleas assembly to the renderer
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(nucleasAssembly);
	renderer->SetBackground(colors->GetColor3d("Black").GetData());

	// Create electron shells and store them into a collection
	vtkSmartPointer<vtkPropCollection> shellCollection = vtkSmartPointer<vtkPropCollection>::New();
	shellCollection = Atom::CreateElectronShell(renderer, shellCollection, 4, colors);

	// Create render window
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetLineSmoothing(true);
	renderWindow->SetPolygonSmoothing(true);
	renderWindow->SetPointSmoothing(true);
	renderWindow->SetMultiSamples(0);
	renderWindow->SetWindowName("Stylized Atom");
	renderWindow->AddRenderer(renderer);

	// Create render window interactor to interact with the actors
	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Switch interactor style to trackball controls
	vtkNew<vtkInteractorStyleSwitch> styleSwitch;
	styleSwitch->SetCurrentStyleToTrackballCamera();
	renderWindowInteractor->SetInteractorStyle(styleSwitch);

	// Create orientation marker to be displayed in the marker widget
#if VTK_HAS_COW
	vtkNew<vtkCameraOrientationWidget> camOrientManipulator;
	camOrientManipulator->SetParentRenderer(renderer);
	// Enable the widget.
	camOrientManipulator->On();
#else
	vtkNew<vtkAxesActor> marker;
	vtkNew<vtkOrientationMarkerWidget> markerWidget;
	double rgba[4]{ 0.0, 0.0, 0.0, 0.0 };
	colors->GetColor("Grey", rgba);
	markerWidget->SetOutlineColor(rgba[0], rgba[1], rgba[2]);
	markerWidget->SetOrientationMarker(marker);
	markerWidget->SetInteractor(renderWindowInteractor);
	markerWidget->SetViewport(0.85, 0.75, 1, 1); // Set the size and position of the widget at the top right corner of the window.
	markerWidget->EnabledOn();
	markerWidget->InteractiveOn();
#endif

	// Set camera position of the renderer in the world coordinates
	vtkCamera* camera = renderer->GetActiveCamera();
	camera->SetPosition(0.0, 0.0, 200.0);
	Atom::currentZoomInFactor = 1;
	renderer->ResetCameraClippingRange();

	// Add legend
	vtkSmartPointer<vtkLegendBoxActor> legendActor = Atom::CreateLegendWidget(colors);
	renderer->AddActor(legendActor);

	// Create an observer to periodically rotate the nucleas
	vtkSmartPointer<vtkCallbackCommand> nucleasCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	nucleasCallback->SetCallback(Atom::RotateNucleas);
	nucleasCallback->SetClientData(nucleasAssembly);
	renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, nucleasCallback);

	// Create an observer to periodically rotate all electron shells
	vtkSmartPointer<vtkCallbackCommand> shellCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	shellCallback->SetCallback(Atom::RotateShell);
	shellCallback->SetClientData(shellCollection);
	renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, shellCallback);

	// Create an observer to adjust the camera zooming at startup
	vtkSmartPointer<vtkCallbackCommand> dollyCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	dollyCallback->SetCallback(Atom::AdjustZoom);
	dollyCallback->SetClientData(camera);
	Atom::dollyObserverTag = renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, dollyCallback);

	//renderWindow->SetFullScreen(true); // Uncomment this line for fullscreen window size
	renderWindow->Render();
	renderWindowInteractor->Initialize();
	renderWindowInteractor->CreateRepeatingTimer(0); // Create repeating timer for the rotation callbacks
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}

namespace Atom {

	vtkSmartPointer<vtkSphereSource> CreateParticleSource(double* center, double radius, double* color) {
		vtkSmartPointer<vtkSphereSource> particleSource = vtkSmartPointer<vtkSphereSource>::New();
		particleSource->SetRadius(radius);
		particleSource->SetPhiResolution(100);
		particleSource->SetThetaResolution(100);
		particleSource->SetCenter(center);
		particleSource->Update();
		
		return particleSource;
	}

	vtkSmartPointer<vtkActor> CreateParticle(double* center, double radius, double* color) {
		vtkSmartPointer<vtkSphereSource> particleSource = CreateParticleSource(center, radius, color);

		vtkSmartPointer<vtkPolyDataMapper> particleMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		particleMapper->SetInputConnection(particleSource->GetOutputPort());

		vtkSmartPointer<vtkActor> particleActor = vtkSmartPointer<vtkActor>::New();
		particleActor->SetMapper(particleMapper);
		particleActor->GetProperty()->SetColor(color);
		return particleActor;
	}

	vtkSmartPointer<vtkAssembly> CreateNucleas(vtkSmartPointer<vtkAssembly> nucleusAssembly,
                                               int nucleasRadius,
                                               int protons,
                                               int neutrons,
                                               vtkSmartPointer<vtkNamedColors> colors) {
		// Randomly position each particle inside a bigger virtual sphere of radius nucleasRadius.
		// A smaller value of nucleusRadius will result in a more spherical nucleus.
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> dist(0, 1);
		std::uniform_real_distribution<double> theta_dist(0.0, 2.0 * vtkMath::Pi());
		std::uniform_real_distribution<double> phi_dist(0.0, vtkMath::Pi());

		int totalParticles = protons + neutrons;
		for (int i = 1; i <= totalParticles; i++) {
			double u = dist(gen);
			double R = nucleasRadius * cbrt(u);

			double theta = theta_dist(gen);
			double phi = phi_dist(gen);
			double center[3]{ R * sin(phi) * cos(theta),R * sin(phi) * sin(theta), R * cos(phi) };

			// Add protons first then neutron
			double* color = (i <= protons) ? colors->GetColor3d("Proton").GetData() : colors->GetColor3d("Neutron").GetData();
			vtkSmartPointer<vtkActor> actor = CreateParticle(center, 1, color);

			vtkProperty* actorProperty = actor->GetProperty();
			actorProperty->SetInterpolationToPBR();
			actorProperty->SetMetallic(0.15);

			nucleusAssembly->AddPart(actor);
		}

		nucleusAssembly->SetObjectName("Nucleas");
		return nucleusAssembly;
	}

	vtkSmartPointer<vtkPropCollection> CreateElectronShell(vtkSmartPointer<vtkRenderer> renderer,
                                                           vtkSmartPointer<vtkPropCollection> shellCollection,
                                                           int shells,
                                                           vtkSmartPointer<vtkNamedColors> colors) {
		double normal[3] = { 0.0, 1.0, 0.3 };
		vtkMath::Normalize(normal);

		double zAxis[3] = { 0.0, 0.0, 1.0 };

		double axisOfRotation[3];
		vtkMath::Cross(zAxis, normal, axisOfRotation);

		double angle = acos(vtkMath::Dot(zAxis, normal)) * 180.0 / vtkMath::Pi(); // Convert radians to degrees

		vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
		transform->RotateWXYZ(angle, axisOfRotation);

		for (int i = 1; i <= shells; i++) {
			vtkSmartPointer<vtkEllipseArcSource> shellSource = vtkSmartPointer<vtkEllipseArcSource>::New();
			shellSource->SetCenter(0.0, 0.0, 0.0);
			shellSource->SetRatio(1.0);
			shellSource->SetMajorRadiusVector(0.0, 0.0, 10.0 + i * 3.0);
			shellSource->SetStartAngle(0);
			shellSource->SetSegmentAngle(360);
			shellSource->SetNormal(normal);
			shellSource->Update();

			vtkSmartPointer<vtkPolyDataMapper> shellMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			shellMapper->SetInputConnection(shellSource->GetOutputPort());

			vtkSmartPointer<vtkActor> shellActor = vtkSmartPointer<vtkActor>::New();
			shellActor->SetMapper(shellMapper);

			vtkProperty* shellProperty = shellActor->GetProperty();
			shellProperty->SetColor(colors->GetColor3d("Shell").GetData());
			shellProperty->SetLineWidth(0.5);

			vtkSmartPointer<vtkAssembly> shellAssembly = vtkSmartPointer<vtkAssembly>::New();
			shellAssembly->AddPart(shellActor);

			// 2n^2 rule
			int totalElectrons = 2 * pow(i, 2);
			for (int j = 0; j < totalElectrons; j++) {
				double angle = 2.0 * vtkMath::Pi() * j / totalElectrons;
				double radius = 10.0 + i * 3;
				double x = radius * cos(angle);
				double y = radius * sin(angle);
				double z = 0;

				double electronPos[3] = { x, y, z };
				transform->TransformPoint(electronPos, electronPos);

				vtkSmartPointer<vtkActor> electronActor = CreateParticle(electronPos, 0.5, colors->GetColor3d("Electron").GetData());
				shellAssembly->AddPart(electronActor);
			}

			renderer->AddActor(shellAssembly);
			shellCollection->AddItem(shellAssembly);
		}

		shellCollection->SetObjectName("ShellCollection");
		return shellCollection;
	}

	void RotateNucleas(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData) {
		vtkRenderWindowInteractor* interactor = static_cast<vtkRenderWindowInteractor*>(caller);
		vtkAssembly* assembly = static_cast<vtkAssembly*>(clientData);

		vtkTransform* transform = vtkTransform::SafeDownCast(assembly->GetUserTransform());
		if (!transform) {
			transform = vtkTransform::New();
			assembly->SetUserTransform(transform);
		}

		transform->RotateY(5.0);
		interactor->GetRenderWindow()->Render();
	}

	void RotateShell(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData) {
		vtkRenderWindowInteractor* interactor = static_cast<vtkRenderWindowInteractor*>(caller);
		vtkPropCollection* assemblyCollection = static_cast<vtkPropCollection*>(clientData);

		assemblyCollection->InitTraversal();
		for (vtkIdType i = 0; i < assemblyCollection->GetNumberOfItems(); i++) {
			vtkAssembly* assembly = vtkAssembly::SafeDownCast(assemblyCollection->GetNextItemAsObject());
			vtkTransform* transform = vtkTransform::SafeDownCast(assembly->GetUserTransform());
			if (!transform) {
				transform = vtkTransform::New();
				assembly->SetUserTransform(transform);
			}

			double rotation[3] = { 0.0, 10.0, 10.0 };
			double angle = (i % 2 != 0) ? 5.0 : -5.0;
			transform->RotateWXYZ(angle, rotation);
			interactor->GetRenderWindow()->Render();
		}
	}

	void AdjustZoom(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData) {
		vtkRenderWindowInteractor* interactor = static_cast<vtkRenderWindowInteractor*>(caller);
		vtkCamera* camera = static_cast<vtkCamera*>(clientData);

		double factor = currentZoomInFactor + ZOOM_IN_DELTA;
		currentZoomInFactor = factor;
		camera->Zoom(factor);

		if (currentZoomInFactor >= 1.16) {
			interactor->RemoveObserver(dollyObserverTag);
		}
		interactor->GetRenderWindow()->Render();
	}

	vtkSmartPointer<vtkLegendBoxActor> CreateLegendWidget(vtkSmartPointer<vtkNamedColors> colors) {
		double center[3] = { 0.0, 0.0, 0.0 };
		auto protonLegend = Atom::CreateParticleSource(center, 1, colors->GetColor3d("Proton").GetData());
		auto neutronLegend = Atom::CreateParticleSource(center, 1, colors->GetColor3d("Neutron").GetData());
		auto electronLegend = Atom::CreateParticleSource(center, 1, colors->GetColor3d("Electron").GetData());
		auto electronShellLegend = vtkSmartPointer<vtkRegularPolygonSource>::New();
		electronShellLegend->GeneratePolygonOff();
		electronShellLegend->SetRadius(1.0);
		electronShellLegend->SetNumberOfSides(100);
		electronShellLegend->SetCenter(center);
		electronShellLegend->Update();

		vtkSmartPointer<vtkLegendBoxActor> legendActor = vtkSmartPointer<vtkLegendBoxActor>::New();
		legendActor->SetHeight(legendActor->GetHeight() - 0.1);
		legendActor->SetWidth(legendActor->GetWidth() - 0.1);

		legendActor->SetNumberOfEntries(4);
		legendActor->SetEntry(2, protonLegend->GetOutput(), "Proton", colors->GetColor3d("Proton").GetData());
		legendActor->SetEntry(0, neutronLegend->GetOutput(), "Neutron", colors->GetColor3d("Neutron").GetData());
		legendActor->SetEntry(1, electronLegend->GetOutput(), "Electron", colors->GetColor3d("Electron").GetData());
		legendActor->SetEntry(3, electronShellLegend->GetOutput(), "Energy Level", colors->GetColor3d("Shell").GetData());

		vtkTextProperty* legendTextProperty = legendActor->GetEntryTextProperty();
		legendTextProperty->SetFontFamilyToTimes();
		//legendTextProperty->SetFrame(true);
		legendTextProperty->SetUseTightBoundingBox(true);
		legendTextProperty->BoldOn();

		legendActor->BorderOff();
		legendActor->GetPositionCoordinate()->SetCoordinateSystemToView();
		legendActor->GetPositionCoordinate()->SetValue(.8, -1.0);
		legendActor->GetPosition2Coordinate()->SetCoordinateSystemToView();
		legendActor->GetPosition2Coordinate()->SetValue(1, -0.7);
		legendActor->UseBackgroundOn();
		legendActor->SetBackgroundColor(colors->GetColor3d("White").GetData());
		
		return legendActor;
	}
}
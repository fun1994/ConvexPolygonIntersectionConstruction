#include "ConvexPolygonIntersectionConstruction.h"
#include "ConvexHull.h"

void print(Point& p) {
	std::cout << "(" << p.x << ", " << p.y << ")";
}

void print(std::vector<Point>& P) {
	for (int i = 0; i < P.size(); i++) {
		print(P[i]);
		std::cout << " ";
	}
	std::cout << std::endl;
}

void generatePoints(std::vector<std::vector<Point>>& S1, std::vector<std::vector<Point>>& S2) {
	Point p;
	S1.push_back(std::vector<Point>());
	p.x = 0.550658222221381; p.y = 0.8730573565870582; S1[0].push_back(p);
	p.x = 0.10019098082769762; p.y = 0.958904283087749; S1[0].push_back(p);
	p.x = 0.05613023961031727; p.y = 0.46915778022759025; S1[0].push_back(p);
	p.x = 0.6198299907311171; p.y = 0.598800802182927; S1[0].push_back(p);
	p.x = 0.5716063392454527; p.y = 0.9261979685131304; S1[0].push_back(p);
	p.x = 0.343847410515538; p.y = 0.38502303779320746; S1[0].push_back(p);
	p.x = 0.524213020085236; p.y = 0.7610381815943781; S1[0].push_back(p);
	p.x = 0.9454455508212812; p.y = 0.7394484167773843; S1[0].push_back(p);
	p.x = 0.7877564181084916; p.y = 0.04148529431967485; S1[0].push_back(p);
	p.x = 0.7163037565183675; p.y = 0.42108638482993244; S1[0].push_back(p);
	S1.push_back(std::vector<Point>());
	p.x = 0.29675635338942197; p.y = 0.8681939956794339; S1[1].push_back(p);
	p.x = 0.28158489257854913; p.y = 0.7615822191269392; S1[1].push_back(p);
	p.x = 0.08351696925926555; p.y = 0.530909753675054; S1[1].push_back(p);
	p.x = 0.8474494118945182; p.y = 0.42159964907504; S1[1].push_back(p);
	p.x = 0.37429865709252064; p.y = 0.7010555950698281; S1[1].push_back(p);
	p.x = 0.3192034758471758; p.y = 0.38907772156328135; S1[1].push_back(p);
	p.x = 0.015064531429036854; p.y = 0.2599019451109652; S1[1].push_back(p);
	p.x = 0.21465899252023624; p.y = 0.07840769958924465; S1[1].push_back(p);
	p.x = 0.027998324460973056; p.y = 0.02161813348403363; S1[1].push_back(p);
	p.x = 0.5501429095592411; p.y = 0.4891937005689839; S1[1].push_back(p);
	S1.push_back(std::vector<Point>());
	p.x = 0.44101435336517203; p.y = 0.44132737958390733; S1[2].push_back(p);
	p.x = 0.742067145226046; p.y = 0.7321839804179926; S1[2].push_back(p);
	p.x = 0.6908565476380871; p.y = 0.7084622177964962; S1[2].push_back(p);
	p.x = 0.8489314760997054; p.y = 0.4324764940354422; S1[2].push_back(p);
	p.x = 0.7058092940979401; p.y = 0.18557762885150064; S1[2].push_back(p);
	p.x = 0.727957862407473; p.y = 0.764182434262949; S1[2].push_back(p);
	p.x = 0.8254621285422171; p.y = 0.43684808351027693; S1[2].push_back(p);
	p.x = 0.8836134989727715; p.y = 0.3117866772852972; S1[2].push_back(p);
	p.x = 0.897296326624907; p.y = 0.3111450125190578; S1[2].push_back(p);
	p.x = 0.6705552362821976; p.y = 0.8257034039676939; S1[2].push_back(p);
	S1.push_back(std::vector<Point>());
	p.x = 0.837197502399392; p.y = 0.6276252407524356; S1[3].push_back(p);
	p.x = 0.48981650062954263; p.y = 0.29644107455060353; S1[3].push_back(p);
	p.x = 0.5563688786276644; p.y = 0.2818858017928104; S1[3].push_back(p);
	p.x = 0.8213629057683982; p.y = 0.3090219893957392; S1[3].push_back(p);
	p.x = 0.8868442052753597; p.y = 0.5614718047665357; S1[3].push_back(p);
	p.x = 0.012792252818687455; p.y = 0.21131528964333357; S1[3].push_back(p);
	p.x = 0.05775082145058419; p.y = 0.47434953984283124; S1[3].push_back(p);
	p.x = 0.3082901818145698; p.y = 0.8438272175893663; S1[3].push_back(p);
	p.x = 0.9385812332583701; p.y = 0.39049399661175466; S1[3].push_back(p);
	p.x = 0.18006143703271482; p.y = 0.3493465613389499; S1[3].push_back(p);
	S1.push_back(std::vector<Point>());
	p.x = 0.14598313429409437; p.y = 0.34164297730768034; S1[4].push_back(p);
	p.x = 0.212047578931297; p.y = 0.757315031941703; S1[4].push_back(p);
	p.x = 0.2544016038057638; p.y = 0.03256379289890887; S1[4].push_back(p);
	p.x = 0.47728234759085564; p.y = 0.3380954953606622; S1[4].push_back(p);
	p.x = 0.2926222287529944; p.y = 0.8487778940966624; S1[4].push_back(p);
	p.x = 0.6258047599324899; p.y = 0.8902705439274289; S1[4].push_back(p);
	p.x = 0.5218664146976408; p.y = 0.5045739213171125; S1[4].push_back(p);
	p.x = 0.7287626439208306; p.y = 0.22891246465750192; S1[4].push_back(p);
	p.x = 0.1521747469156174; p.y = 0.9489995045272419; S1[4].push_back(p);
	p.x = 0.7125818318131212; p.y = 0.8205958896722562; S1[4].push_back(p);
	S1.push_back(std::vector<Point>());
	p.x = 0.9954743365587353; p.y = 0.6353746913847197; S1[5].push_back(p);
	p.x = 0.48038535522670667; p.y = 0.2478836752570699; S1[5].push_back(p);
	p.x = 0.8630683925743325; p.y = 0.3020359107888325; S1[5].push_back(p);
	p.x = 0.3780206656810694; p.y = 0.06948308102859546; S1[5].push_back(p);
	p.x = 0.7025248700734514; p.y = 0.79721370077534; S1[5].push_back(p);
	p.x = 0.22979501832595262; p.y = 0.8373522889223982; S1[5].push_back(p);
	p.x = 0.7821070460550387; p.y = 0.3597048483944624; S1[5].push_back(p);
	p.x = 0.728603646862922; p.y = 0.28050977872434024; S1[5].push_back(p);
	p.x = 0.8376389865096059; p.y = 0.11509260581246183; S1[5].push_back(p);
	p.x = 0.7348010107813561; p.y = 0.995083212042214; S1[5].push_back(p);
	S1.push_back(std::vector<Point>());
	p.x = 0.44812825784456334; p.y = 0.4230369963547266; S1[6].push_back(p);
	p.x = 0.2805226517021949; p.y = 0.1053165575170123; S1[6].push_back(p);
	p.x = 0.9306095184674695; p.y = 0.4081265494258052; S1[6].push_back(p);
	p.x = 0.9640842748803236; p.y = 0.644852589635744; S1[6].push_back(p);
	p.x = 0.9672058394155697; p.y = 0.38105094703306075; S1[6].push_back(p);
	p.x = 0.36670985162647485; p.y = 0.4653728437623651; S1[6].push_back(p);
	p.x = 0.5789174920814465; p.y = 0.5758208450892132; S1[6].push_back(p);
	p.x = 0.4235397773506594; p.y = 0.031427612284111506; S1[6].push_back(p);
	p.x = 0.9756228409984515; p.y = 0.42386793378931054; S1[6].push_back(p);
	p.x = 0.11094773945194347; p.y = 0.7294156104202699; S1[6].push_back(p);
	S1.push_back(std::vector<Point>());
	p.x = 0.6954239633778929; p.y = 0.898281550318704; S1[7].push_back(p);
	p.x = 0.12015573750221953; p.y = 0.9138527657846174; S1[7].push_back(p);
	p.x = 0.8410014239045521; p.y = 0.18548494429140716; S1[7].push_back(p);
	p.x = 0.5687119519792214; p.y = 0.9914903397670217; S1[7].push_back(p);
	p.x = 0.8400673702720415; p.y = 0.9673662861581611; S1[7].push_back(p);
	p.x = 0.9185452728826294; p.y = 0.8157285559899213; S1[7].push_back(p);
	p.x = 0.25780419152588807; p.y = 0.3378128489219193; S1[7].push_back(p);
	p.x = 0.5133435549198113; p.y = 0.08836441338580414; S1[7].push_back(p);
	p.x = 0.6609614456352272; p.y = 0.419026339693097; S1[7].push_back(p);
	p.x = 0.5150256216092122; p.y = 0.9069554104310675; S1[7].push_back(p);
	S1.push_back(std::vector<Point>());
	p.x = 0.4374661537466896; p.y = 0.4340583330179689; S1[8].push_back(p);
	p.x = 0.052867167339262555; p.y = 0.7906040010596511; S1[8].push_back(p);
	p.x = 0.41203371510591646; p.y = 0.7299113989896239; S1[8].push_back(p);
	p.x = 0.9551168559977239; p.y = 0.6705652598148538; S1[8].push_back(p);
	p.x = 0.4952969646199411; p.y = 0.012521272718835408; S1[8].push_back(p);
	p.x = 0.9490811425452127; p.y = 0.7152114066235483; S1[8].push_back(p);
	p.x = 0.5762376711897779; p.y = 0.360502506643565; S1[8].push_back(p);
	p.x = 0.9164077499294796; p.y = 0.8354406404168839; S1[8].push_back(p);
	p.x = 0.8879770498154577; p.y = 0.7678408068126006; S1[8].push_back(p);
	p.x = 0.7994703063014785; p.y = 0.19303686650964202; S1[8].push_back(p);
	S1.push_back(std::vector<Point>());
	p.x = 0.005760338491872852; p.y = 0.509341923908721; S1[9].push_back(p);
	p.x = 0.23161268835022097; p.y = 0.9985398882203657; S1[9].push_back(p);
	p.x = 0.9346538613720049; p.y = 0.972559417634858; S1[9].push_back(p);
	p.x = 0.4247675922100522; p.y = 0.05361261728328104; S1[9].push_back(p);
	p.x = 0.5625838247674789; p.y = 0.48581017565889195; S1[9].push_back(p);
	p.x = 0.36858218930131037; p.y = 0.7333250308928659; S1[9].push_back(p);
	p.x = 0.7807774028909151; p.y = 0.1541286260332304; S1[9].push_back(p);
	p.x = 0.26867934758049494; p.y = 0.2690282042802411; S1[9].push_back(p);
	p.x = 0.42256280428193793; p.y = 0.23825573420107216; S1[9].push_back(p);
	p.x = 0.3101256758048948; p.y = 0.2492227224889425; S1[9].push_back(p);
	S2.push_back(std::vector<Point>());
	p.x = 0.4035363780955925; p.y = 0.010911752848657974; S2[0].push_back(p);
	p.x = 0.7670163042936525; p.y = 0.5041836463207697; S2[0].push_back(p);
	p.x = 0.3843681626112472; p.y = 0.8880247928610145; S2[0].push_back(p);
	p.x = 0.06811604817434658; p.y = 0.17902318523494565; S2[0].push_back(p);
	p.x = 0.17956182075001725; p.y = 0.2580174883028983; S2[0].push_back(p);
	p.x = 0.8813635533289664; p.y = 0.025725294458332826; S2[0].push_back(p);
	p.x = 0.5970903737306514; p.y = 0.16604665891865167; S2[0].push_back(p);
	p.x = 0.07041797849819653; p.y = 0.27471472108497663; S2[0].push_back(p);
	p.x = 0.3504882661702797; p.y = 0.7174291510136083; S2[0].push_back(p);
	p.x = 0.23157361613974292; p.y = 0.41229085214889516; S2[0].push_back(p);
	S2.push_back(std::vector<Point>());
	p.x = 0.5386856939784389; p.y = 0.2582980394499207; S2[1].push_back(p);
	p.x = 0.9805005938123325; p.y = 0.5269937852967587; S2[1].push_back(p);
	p.x = 0.20784957257641234; p.y = 0.7969368207232156; S2[1].push_back(p);
	p.x = 0.7334835296600671; p.y = 0.383868823561762; S2[1].push_back(p);
	p.x = 0.30400213897085093; p.y = 0.0980378614994376; S2[1].push_back(p);
	p.x = 0.41370491748083627; p.y = 0.18459020695768102; S2[1].push_back(p);
	p.x = 0.4449920212198899; p.y = 0.16108351663503062; S2[1].push_back(p);
	p.x = 0.2239172082347496; p.y = 0.7613317037304301; S2[1].push_back(p);
	p.x = 0.15580683161553832; p.y = 0.35347821254930345; S2[1].push_back(p);
	p.x = 0.22783364550128182; p.y = 0.2704280848853101; S2[1].push_back(p);
	S2.push_back(std::vector<Point>());
	p.x = 0.16270710651499043; p.y = 0.36960260547275725; S2[2].push_back(p);
	p.x = 0.8344665654383839; p.y = 0.9884271200991699; S2[2].push_back(p);
	p.x = 0.44012945713413354; p.y = 0.9738114984387817; S2[2].push_back(p);
	p.x = 0.7626954350252338; p.y = 0.0980086132004987; S2[2].push_back(p);
	p.x = 0.2976932855433174; p.y = 0.19989674773325083; S2[2].push_back(p);
	p.x = 0.10559080607019211; p.y = 0.8197714535504762; S2[2].push_back(p);
	p.x = 0.27123198071910104; p.y = 0.043116000221186335; S2[2].push_back(p);
	p.x = 0.07252456360469273; p.y = 0.7573203391863603; S2[2].push_back(p);
	p.x = 0.8342150616534556; p.y = 0.015677601606404967; S2[2].push_back(p);
	p.x = 0.7564881764557724; p.y = 0.9211322228266203; S2[2].push_back(p);
	S2.push_back(std::vector<Point>());
	p.x = 0.7902520073180519; p.y = 0.8788313340904775; S2[3].push_back(p);
	p.x = 0.4539067766407966; p.y = 0.018729554754583377; S2[3].push_back(p);
	p.x = 0.812419362669611; p.y = 0.11035393282232075; S2[3].push_back(p);
	p.x = 0.37262363033362067; p.y = 0.9327600263680171; S2[3].push_back(p);
	p.x = 0.5224075623951298; p.y = 0.445193852046295; S2[3].push_back(p);
	p.x = 0.45658444935360454; p.y = 0.12200587556349973; S2[3].push_back(p);
	p.x = 0.2954391403455918; p.y = 0.4709661051897349; S2[3].push_back(p);
	p.x = 0.8172205395357743; p.y = 0.30732323461967537; S2[3].push_back(p);
	p.x = 0.39986499899589434; p.y = 0.04706730284507399; S2[3].push_back(p);
	p.x = 0.09764470098041267; p.y = 0.48155580942562104; S2[3].push_back(p);
	S2.push_back(std::vector<Point>());
	p.x = 0.6598988286471419; p.y = 0.11417494078143176; S2[4].push_back(p);
	p.x = 0.41130306345578616; p.y = 0.32564476133337916; S2[4].push_back(p);
	p.x = 0.7996984032019707; p.y = 0.5706639244529541; S2[4].push_back(p);
	p.x = 0.22034143934096528; p.y = 0.40902534155645975; S2[4].push_back(p);
	p.x = 0.05287837888827718; p.y = 0.5310174527831294; S2[4].push_back(p);
	p.x = 0.5553183834069239; p.y = 0.5341670271615333; S2[4].push_back(p);
	p.x = 0.04373964535713404; p.y = 0.04739767077379142; S2[4].push_back(p);
	p.x = 0.30863260049760755; p.y = 0.10727376670007915; S2[4].push_back(p);
	p.x = 0.09139490907765835; p.y = 0.30051992543898864; S2[4].push_back(p);
	p.x = 0.08643537087945452; p.y = 0.5864428150247146; S2[4].push_back(p);
	S2.push_back(std::vector<Point>());
	p.x = 0.8694658628597248; p.y = 0.11036941147641588; S2[5].push_back(p);
	p.x = 0.727073133034394; p.y = 0.879750607159667; S2[5].push_back(p);
	p.x = 0.7845673502307552; p.y = 0.26735496921400015; S2[5].push_back(p);
	p.x = 0.22788944586896576; p.y = 0.04417167853165049; S2[5].push_back(p);
	p.x = 0.7973113690795457; p.y = 0.9962372777232934; S2[5].push_back(p);
	p.x = 0.7364305085504157; p.y = 0.7144979296246058; S2[5].push_back(p);
	p.x = 0.2442775851134673; p.y = 0.8149376762639429; S2[5].push_back(p);
	p.x = 0.16429840653401573; p.y = 0.9654360163120385; S2[5].push_back(p);
	p.x = 0.06847989521639997; p.y = 0.5748339696720444; S2[5].push_back(p);
	p.x = 0.7924823350386573; p.y = 0.9974087428853938; S2[5].push_back(p);
	S2.push_back(std::vector<Point>());
	p.x = 0.09417636402596119; p.y = 0.3684112346851758; S2[6].push_back(p);
	p.x = 0.19489740915907083; p.y = 0.41888154688490253; S2[6].push_back(p);
	p.x = 0.12024815269322464; p.y = 0.3966357701631744; S2[6].push_back(p);
	p.x = 0.55940046136494; p.y = 0.8657899107470672; S2[6].push_back(p);
	p.x = 0.08958528489040518; p.y = 0.16072632262669628; S2[6].push_back(p);
	p.x = 0.5758849606207938; p.y = 0.4594345690461745; S2[6].push_back(p);
	p.x = 0.13699495793014782; p.y = 0.2986084183688642; S2[6].push_back(p);
	p.x = 0.36169919372209536; p.y = 0.5761810141619431; S2[6].push_back(p);
	p.x = 0.2218005177601674; p.y = 0.15392000025748698; S2[6].push_back(p);
	p.x = 0.6716162647580709; p.y = 0.6104192802518419; S2[6].push_back(p);
	S2.push_back(std::vector<Point>());
	p.x = 0.15786339309616404; p.y = 0.744812475843163; S2[7].push_back(p);
	p.x = 0.114182359861696; p.y = 0.7386574042484958; S2[7].push_back(p);
	p.x = 0.4438749915143253; p.y = 0.27235484747887506; S2[7].push_back(p);
	p.x = 0.6878725498706488; p.y = 0.7850123468580091; S2[7].push_back(p);
	p.x = 0.03189850229445801; p.y = 0.8149473640470593; S2[7].push_back(p);
	p.x = 0.01918882532988264; p.y = 0.13510435267847953; S2[7].push_back(p);
	p.x = 0.14099605265175585; p.y = 0.4098105254419374; S2[7].push_back(p);
	p.x = 0.8688353942144578; p.y = 0.1550721914186033; S2[7].push_back(p);
	p.x = 0.7040444648883173; p.y = 0.5442588073274021; S2[7].push_back(p);
	p.x = 0.7452709454817278; p.y = 0.27204758477514823; S2[7].push_back(p);
	S2.push_back(std::vector<Point>());
	p.x = 0.8701021226177909; p.y = 0.5195259862990513; S2[8].push_back(p);
	p.x = 0.6757271982953821; p.y = 0.7257104613362187; S2[8].push_back(p);
	p.x = 0.4977737180880947; p.y = 0.32370156460835264; S2[8].push_back(p);
	p.x = 0.10180084863473182; p.y = 0.8129682052834798; S2[8].push_back(p);
	p.x = 0.5792065371663723; p.y = 0.7932989024712552; S2[8].push_back(p);
	p.x = 0.22121992305542826; p.y = 0.9256256007253112; S2[8].push_back(p);
	p.x = 0.9962506206923292; p.y = 0.45387994693995104; S2[8].push_back(p);
	p.x = 0.9977699588195024; p.y = 0.6274676138549876; S2[8].push_back(p);
	p.x = 0.3554081298370876; p.y = 0.48742556599745035; S2[8].push_back(p);
	p.x = 0.5345930727370231; p.y = 0.7204341986814766; S2[8].push_back(p);
	S2.push_back(std::vector<Point>());
	p.x = 0.5545813927359364; p.y = 0.2856695785336564; S2[9].push_back(p);
	p.x = 0.15317806301876946; p.y = 0.8259961821084734; S2[9].push_back(p);
	p.x = 0.6423931786930728; p.y = 0.5724621272166783; S2[9].push_back(p);
	p.x = 0.5427427114343433; p.y = 0.8619496034096388; S2[9].push_back(p);
	p.x = 0.3705191332326665; p.y = 0.8902037158913209; S2[9].push_back(p);
	p.x = 0.06885460093032714; p.y = 0.6568711996642632; S2[9].push_back(p);
	p.x = 0.5213177095209287; p.y = 0.4129491577014527; S2[9].push_back(p);
	p.x = 0.5530652707479945; p.y = 0.2406795513166584; S2[9].push_back(p);
	p.x = 0.9367185664196244; p.y = 0.28510606583394993; S2[9].push_back(p);
	p.x = 0.3444998909204149; p.y = 0.8502903204525539; S2[9].push_back(p);
}

bool isSame(std::vector<std::vector<Point>>& P) {
	for (int i = 1; i < P.size(); i++) {
		if (P[0].size() != P[i].size()) {
			return false;
		}
		int index = -1;
		for (int j = 0; j < P[0].size(); j++) {
			if (abs(P[0][0].x - P[i][j].x) < 1e-12 && abs(P[0][0].y - P[i][j].y) < 1e-12) {
				index = j;
				break;
			}
		}
		if (!P[0].empty() && index == -1) {
			return false;
		}
		for (int j = 0; j < P[0].size(); j++) {
			if (abs(P[0][j].x - P[i][(index + j) % P[0].size()].x) > 1e-12 || abs(P[0][j].y - P[i][(index + j) % P[0].size()].y) > 1e-12) {
				return false;
			}
		}
	}
	return true;
}

void test() {
	ConvexPolygonIntersectionConstruction CPIC;
	ConvexHull CH;
	std::vector<std::vector<Point>> S1, S2;
	generatePoints(S1, S2);
	for (int i = 0; i < S1.size(); i++) {
		std::vector<Point> P1 = CH.GrahamScan(S1[i]);
		std::vector<Point> P2 = CH.GrahamScan(S2[i]);
		for (int j = 0; j < P1.size(); j++) {
			P1.push_back(P1[0]);
			P1.erase(P1.begin());
			for (int k = 0; k < P2.size(); k++) {
				P2.push_back(P2[0]);
				P2.erase(P2.begin());
				for (int m = 0; m < 20; m++) {
					for (int p = 0; p < P2.size(); p++) {
						P2[p].x += 0.1;
					}
					for (int n = 0; n < 20; n++) {
						for (int p = 0; p < P2.size(); p++) {
							P2[p].y += 0.1;
						}
						std::vector<std::vector<Point>> P;
						P.push_back(CPIC.bruteForce(P1, P2));
						P.push_back(CPIC.ORourke(P1, P2));
						P.push_back(CPIC.planeSweeping(P1, P2));
						if (!isSame(P)) {
							print(P1);
							print(P2);
							std::cin.get();
						}
					}
					for (int p = 0; p < P2.size(); p++) {
						P2[p].y -= 2;
					}
				}
				for (int p = 0; p < P2.size(); p++) {
					P2[p].x -= 2;
				}
			}
		}
	}
}

int main() {
	test();
	return 0;
}

function interpolated = RhieChow( U_cellVar,V_cellVar,U_aP,V_aP,P,DPx,DPy,U_BC,V_BC )
%RHIECHOW Summary of this function goes here
%   Detailed explanation goes here

N = P.domain.dims;

%U_aP_sum = reshape(U_aP_sum, N+2);
%V_aP_sum = reshape(V_aP_sum, N+2);
%U_aP_sumCellVar.value=U_aP_sum;
%U_aP_sumCellVar.domain=P.domain;
%V_aP_sumCellVar.value=V_aP_sum;
%V_aP_sumCellVar.domain=P.domain;

U_cellVarArithmeticFaceValue=arithmeticMean(U_cellVar);
V_cellVarArithmeticFaceValue=arithmeticMean(V_cellVar);

apU_ArithmeticFaceValue=arithmeticMean(U_aP);
apV_ArithmeticFaceValue=arithmeticMean(V_aP);

DPx_ArithmeticFaceValue=arithmeticMean(DPx);
DPx_ArithmeticFaceValue.xvalue(:,[1 end])=0;
DPy_ArithmeticFaceValue=arithmeticMean(DPy);
DPy_ArithmeticFaceValue.yvalue([1 end],:)=0;

facepresgrad=gradientTerm(P);
facepresgrad.xvalue(:,[1 end])=0;
facepresgrad.yvalue([1 end],:)=0;

% from http://www.ctcms.nist.gov/fipy/examples/flow/generated/examples.flow.stokesCavity.html
ue1=U_cellVarArithmeticFaceValue.xvalue+P.domain.cellsize.x(1)*P.domain.cellsize.y(1)./apU_ArithmeticFaceValue.xvalue.*(DPx_ArithmeticFaceValue.xvalue-facepresgrad.xvalue);
% set BC
%ue1
un1=V_cellVarArithmeticFaceValue.yvalue+P.domain.cellsize.x(1)*P.domain.cellsize.y(1)./apV_ArithmeticFaceValue.yvalue.*(DPy_ArithmeticFaceValue.yvalue-facepresgrad.yvalue);

ue1(1,:)=0;
ue1(end,:)=0;
%ue1(2:end-1,1)=1;

un1(:,1)=0;
un1(:,end)=0;

interpolated.xvalue=ue1;
interpolated.yvalue=un1;
interpolated.domain=U_cellVar.domain;


%test=createFaceVariable(P.domain,ue1,U_BC);


%ue=(U_cellVar.value(1:end-1,:)+U_cellVar.value(2:end,:))/2-0.5*(1./U_aP_sum(1:end-1,:)+1./U_aP_sum(2:end,:)).*(P.value(2:end,:)-P.value(1:end-1,:))+0.5*(DPx.value(1:end-1,:)./U_aP_sum(1:end-1,:)+DPx.value(2:end,:)./U_aP_sum(2:end,:));
%ue=(U_cellVar.value(1:end-1,:)+U_cellVar.value(2:end,:))/2-0.5*P.domain.cellsize.y(1)*(1./U_aP_sum(1:end-1,:)+1./U_aP_sum(2:end,:)).*(P.value(2:end,:)-P.value(1:end-1,:))-0.5./U_aP_sum.*(DPx.value(1:end-1,:)+DPx.value(2:end,:));

%un=(V_cellVar.value(:,1:end-1)+V_cellVar.value(:,2:end))/2-0.5*(1./U_aP_sum(:,1:end-1)+1./U_aP_sum(:,2:end)).*(P.value(:,2:end)-P.value(:,1:end-1))+0.5*(DPy.value(:,1:end-1)./U_aP_sum(:,1:end-1)+DPy.value(:,2:end)./U_aP_sum(:,2:end));
%un=(V_cellVar.value(:,1:end-1)+V_cellVar.value(:,2:end))/2-0.5*P.domain.cellsize.x(1)*(1./U_aP_sum(:,1:end-1)+1./U_aP_sum(:,2:end)).*(P.value(:,2:end)-P.value(:,1:end-1))+0.5/V_aP_sum.*(DPy.value(:,1:end-1)+DPy.value(:,2:end));

%interpolated.xvalue=ue(:,2:end-1);
%interpolated.yvalue=un(2:end-1,:);
%interpolated.zvalue=[];
%interpolated.domain=U_cellVar.domain;

end

